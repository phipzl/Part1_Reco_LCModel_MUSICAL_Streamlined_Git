#!/usr/bin/env python3
"""Warn when the T1 given with -t and the MRSI are in different frames of reference.

    check_t1_frame.py <T1 DICOM folder or NIfTI/MINC file> <MRSI twix .dat>

The T1-based mask is placed through scanner coordinates, which only match within one frame of
reference: a new landmark between the two scans shifts the mask. Warns only, always exits 0.
"""

import datetime
import json
import os
import re
import struct
import sys

ADVICE = (
    "  The T1-based mask, segmentation and overlays are shifted relative to the MRSI.\n"
    "  Use a T1 acquired after the last landmark, or register the T1 to the MRSI before\n"
    "  passing it with -t."
)
MAX_HEADER_BYTES = 64 * 1024**2

FRAME_OF_REFERENCE_UID = (0x0020, 0x0052)
TRANSFER_SYNTAX_UID = (0x0002, 0x0010)
ITEM_DELIMITER = (0xFFFE, 0xE00D)
SEQUENCE_DELIMITER = (0xFFFE, 0xE0DD)
UNDEFINED_LENGTH = 0xFFFFFFFF
IMPLICIT_LITTLE_ENDIAN = "1.2.840.10008.1.2"
UNREADABLE_SYNTAXES = ("1.2.840.10008.1.2.1.99", "1.2.840.10008.1.2.2")  # deflated, big endian
# Explicit VRs with a 2-byte length; all others have 2 reserved bytes and a 4-byte length
SHORT_VRS = b"AE AS AT CS DA DS DT FL FD IS LO LT PN SH SL SS ST TM UI UL US".split()


def twix_frame_of_reference(path):
    """The frame of reference UID in the header of a twix file's last measurement."""
    with open(path, "rb") as f:
        first, n_meas = struct.unpack("<2I", f.read(8))
        offset = 0
        if first == 0 and 0 < n_meas <= 64:  # VD and later: a table of measurements
            f.seek(8 + (n_meas - 1) * 152 + 8)
            (offset,) = struct.unpack("<Q", f.read(8))
        f.seek(offset)
        (header_length,) = struct.unpack("<I", f.read(4))
        if not 0 < header_length <= MAX_HEADER_BYTES:
            raise ValueError("the MRSI file has no twix header")
        header = f.read(header_length - 4).decode("latin1")
    m = re.search(r'"tFrameOfReference">\s*\{[^}]*?"([0-9.]+)"\s*\}', header)
    if m is None:
        raise ValueError("the MRSI header has no frame of reference")
    return m.group(1)


def _element_header(f, implicit):
    """Tag, VR and value length of the next DICOM data element; None at the end of the file."""
    head = f.read(8)
    if len(head) < 8:
        return None
    tag = struct.unpack("<2H", head[:4])
    vr = head[4:6]
    if tag[0] == 0xFFFE or (implicit and tag[0] != 0x0002):
        return tag, b"", struct.unpack("<I", head[4:])[0]
    if vr in SHORT_VRS:
        return tag, vr, struct.unpack("<H", head[6:])[0]
    return tag, vr, struct.unpack("<I", f.read(4))[0]


def _skip_undefined_length(f, implicit, end=SEQUENCE_DELIMITER):
    """Move past a sequence or item of undefined length, up to and including its delimiter."""
    while True:
        element = _element_header(f, implicit)
        if element is None or element[0] == end:
            return
        if element[2] != UNDEFINED_LENGTH:
            f.seek(element[2], 1)
        elif end == SEQUENCE_DELIMITER:  # an item
            _skip_undefined_length(f, implicit, ITEM_DELIMITER)
        else:  # a nested sequence
            _skip_undefined_length(f, implicit)


def dicom_frame_of_reference(path):
    """FrameOfReferenceUID (0020,0052) of a DICOM file, None when it has none."""
    with open(path, "rb") as f:
        f.seek(128)
        if f.read(4) != b"DICM":
            raise ValueError("not a DICOM file")
        implicit = False
        while True:
            element = _element_header(f, implicit)
            if element is None or element[0] > FRAME_OF_REFERENCE_UID:
                return None
            tag, _, length = element
            if length == UNDEFINED_LENGTH:
                _skip_undefined_length(f, implicit)
            elif tag in (TRANSFER_SYNTAX_UID, FRAME_OF_REFERENCE_UID):
                value = f.read(length).decode("latin1").strip(" \0")
                if tag == FRAME_OF_REFERENCE_UID:
                    return value
                if value in UNREADABLE_SYNTAXES:
                    raise ValueError("DICOM transfer syntax " + value + " is not read")
                implicit = value == IMPLICIT_LITTLE_ENDIAN
            else:
                f.seek(length, 1)


def folder_frame_of_reference(folder):
    """FrameOfReferenceUID of the first DICOM file below a folder."""
    for root, dirs, files in os.walk(folder):
        dirs.sort()
        for name in sorted(files):
            try:
                uid = dicom_frame_of_reference(os.path.join(root, name))
            except (OSError, ValueError, struct.error):
                continue
            if uid:
                return uid
            raise ValueError("the T1 DICOM has no FrameOfReferenceUID")
    raise ValueError("no DICOM file in the T1 folder")


def uid_time(uid):
    """Creation time (to the ms) that a Siemens frame of reference UID carries, None without."""
    for part in uid.split("."):
        if len(part) == 17 and part.isdigit():
            try:
                return datetime.datetime.strptime(part, "%Y%m%d%H%M%S%f")
            except ValueError:
                pass
    return None


def same_frame(uid_a, uid_b):
    """Equal UIDs, or equal creation times: the raw data UID differs from the image UID in the
    field before the time even within one frame."""
    return uid_a == uid_b or (uid_time(uid_a) is not None and uid_time(uid_a) == uid_time(uid_b))


def sidecar_time(nifti_path, frame_created):
    """Acquisition time in the JSON sidecar of a NIfTI, on the frame's day if it has no date."""
    sidecar = re.sub(r"\.nii(\.gz)?$", "", nifti_path) + ".json"
    try:
        with open(sidecar) as f:
            meta = json.load(f)
        time = datetime.time.fromisoformat(str(meta["AcquisitionTime"]).split(".")[0])
        day = frame_created.date()
        if "AcquisitionDateTime" in meta:
            day = datetime.date.fromisoformat(str(meta["AcquisitionDateTime"])[:10])
        return datetime.datetime.combine(day, time)
    except (OSError, ValueError, KeyError):
        return None


def frame_report(t1_path, csi_path):
    """The text to print for this T1 and MRSI file, empty when the frames agree."""
    if not os.path.isfile(csi_path):  # -c may hold several files
        csi_path = (csi_path.split() or [""])[0]
    if not csi_path.lower().endswith(".dat"):
        raise ValueError("the MRSI file is not a twix .dat")
    mrsi_uid = twix_frame_of_reference(csi_path)
    mrsi_created = uid_time(mrsi_uid)

    if os.path.isdir(t1_path):
        t1_uid = folder_frame_of_reference(t1_path)
        if same_frame(t1_uid, mrsi_uid):
            return ""
        t1_created = uid_time(t1_uid)
        created = (
            "  T1:   {} (created {:%Y-%m-%d %H:%M:%S})\n  MRSI: {} (created {:%Y-%m-%d %H:%M:%S})\n"
        )
        if t1_created is None or mrsi_created is None:
            created = "  T1:   {0}\n  MRSI: {2}\n"
        return (
            "\nWARNING: the T1 and the MRSI are in different frames of reference:\n"
            + created.format(t1_uid, t1_created, mrsi_uid, mrsi_created)
            + "  The patient was probably re-landmarked in between.\n"
            + ADVICE
            + "\n"
        )

    t1_acquired = mrsi_created and sidecar_time(t1_path, mrsi_created)
    if t1_acquired and t1_acquired < mrsi_created:
        clock = "{:%H:%M:%S}"
        if t1_acquired.date() != mrsi_created.date():
            clock = "{:%Y-%m-%d %H:%M:%S}"
        return (
            "\nWARNING: T1 acquired at "
            + clock.format(t1_acquired)
            + ", before the MRSI frame of reference was created at "
            + clock.format(mrsi_created)
            + ":\n  the patient was probably re-landmarked in between.\n"
            + ADVICE
            + "\n"
        )
    return (
        "T1 frame of reference not compared: a NIfTI or MINC T1 carries no frame of reference UID"
        " (give the T1 DICOM folder with -t to have it compared)"
    )


def main(argv):
    try:
        if len(argv) != 3:
            raise ValueError("usage: check_t1_frame.py <T1> <MRSI .dat>")
        report = frame_report(argv[1], argv[2])
    except Exception as e:  # the check must never stop the pipeline
        report = "T1 frame of reference not checked: {}".format(e)
    if report:
        print(report, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
