#!/usr/bin/env python3
"""Tests of check_t1_frame.py, on synthetic headers only.

python -m unittest discover -s tests
"""

import contextlib
import io
import json
import os
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import check_t1_frame as check  # noqa: E402

EARLY = "1.3.12.2.1107.5.2.0.99999.2.20200102101500123.0.0.0"
LATE = "1.3.12.2.1107.5.2.0.99999.1.20200102110000456.0.0.0"


def twix_header(uid, visible=False):
    text = '<ParamLong."lEchoSpacing">{ }\n<ParamString."tFrameOfReference">{ %s"%s" }\n'
    text = text % ('<Visible> "true" ' if visible else "", uid)
    return struct.pack("<I", len(text) + 4) + text.encode()


def write_twix(path, uid, n_meas=1, **kwargs):
    """A VD-style file: measurement table, then one header per measurement."""
    headers = [twix_header("9.9.9")] * (n_meas - 1) + [twix_header(uid, **kwargs)]
    table, offset = struct.pack("<2I", 0, n_meas), 10240
    for header in headers:
        table += struct.pack("<2I2Q", 1, 1, offset, len(header)) + bytes(128)
        offset += 512
    with open(path, "wb") as f:
        f.write(table)
        for i, header in enumerate(headers):
            f.seek(10240 + 512 * i)
            f.write(header)


def element(group, number, vr, value, implicit=False):
    value += b"\0" * (len(value) % 2)
    tag = struct.pack("<2H", group, number)
    if implicit:
        return tag + struct.pack("<I", len(value)) + value
    if vr in check.SHORT_VRS:
        return tag + vr + struct.pack("<H", len(value)) + value
    return tag + vr + bytes(2) + struct.pack("<I", len(value)) + value


def write_dicom(path, uid, implicit=False):
    """Meta header, a sequence of undefined length holding a decoy UID, then the UID."""
    syntax = b"1.2.840.10008.1.2" if implicit else b"1.2.840.10008.1.2.1"
    undefined = struct.pack("<I", 0xFFFFFFFF)
    decoy = element(0x0020, 0x0052, b"UI", b"9.9.9", implicit)
    item = struct.pack("<2H", 0xFFFE, 0xE000) + undefined + decoy
    item += struct.pack("<2HI", 0xFFFE, 0xE00D, 0)
    sequence = struct.pack("<2H", 0x0008, 0x1140) + (b"" if implicit else b"SQ" + bytes(2))
    sequence += undefined + item + struct.pack("<2HI", 0xFFFE, 0xE0DD, 0)
    body = element(0x0008, 0x0060, b"CS", b"MR", implicit) + sequence
    if uid:
        body += element(0x0020, 0x0052, b"UI", uid.encode(), implicit)
    body += element(0x0028, 0x0010, b"US", struct.pack("<H", 4), implicit)
    with open(path, "wb") as f:
        f.write(bytes(128) + b"DICM" + element(0x0002, 0x0010, b"UI", syntax) + body)


class FrameCheck(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.dat = self.path("mrsi.dat")
        self.t1 = self.path("t1")
        os.mkdir(self.t1)

    def path(self, name):
        return os.path.join(self.dir.name, name)

    def nifti(self, **sidecar):
        if sidecar:
            with open(self.path("t1.json"), "w") as f:
                json.dump(sidecar, f)
        return self.path("t1.nii.gz")

    def test_twix_reads_the_last_measurement(self):
        write_twix(self.dat, LATE, n_meas=2)
        self.assertEqual(check.twix_frame_of_reference(self.dat), LATE)

    def test_twix_value_after_an_attribute(self):
        write_twix(self.dat, LATE, visible=True)
        self.assertEqual(check.twix_frame_of_reference(self.dat), LATE)

    def test_twix_without_measurement_table(self):
        with open(self.dat, "wb") as f:
            f.write(twix_header(LATE))
        self.assertEqual(check.twix_frame_of_reference(self.dat), LATE)

    def test_dicom_uid_explicit_and_implicit(self):
        for implicit in (False, True):
            write_dicom(self.path("a.dcm"), EARLY, implicit)
            self.assertEqual(check.dicom_frame_of_reference(self.path("a.dcm")), EARLY)

    def test_uid_time(self):
        self.assertEqual(str(check.uid_time(EARLY)), "2020-01-02 10:15:00.123000")
        self.assertIsNone(check.uid_time("1.2.3"))

    def test_same_frame_is_silent(self):
        write_twix(self.dat, LATE)
        # an image UID differs from the raw data UID in the field before the time
        for uid in (LATE, LATE.replace(".1.2020", ".2.2020")):
            write_dicom(os.path.join(self.t1, "1.dcm"), uid)
            self.assertEqual(check.frame_report(self.t1, self.dat), "")

    def test_other_frame_without_times_warns(self):
        write_twix(self.dat, "1.2.3")
        write_dicom(os.path.join(self.t1, "1.dcm"), "1.2.4")
        report = check.frame_report(self.t1, self.dat)
        self.assertIn("  T1:   1.2.4\n  MRSI: 1.2.3\n", report)

    def test_other_frame_warns_with_both_uids(self):
        write_twix(self.dat, LATE)
        with open(os.path.join(self.t1, "0_notes.txt"), "w") as f:
            f.write("not a DICOM file")
        write_dicom(os.path.join(self.t1, "1.dcm"), EARLY)
        report = check.frame_report(self.t1, self.dat)
        self.assertTrue(report.startswith("\nWARNING:"))
        for text in (EARLY, LATE, "10:15:00", "11:00:00", "register the T1"):
            self.assertIn(text, report)

    def test_first_of_several_csi_files_is_read(self):
        write_twix(self.dat, LATE)
        write_dicom(os.path.join(self.t1, "1.dcm"), LATE)
        self.assertEqual(check.frame_report(self.t1, self.dat + " " + self.path("b.dat")), "")

    def test_nifti_acquired_before_the_frame_warns(self):
        write_twix(self.dat, LATE)
        report = check.frame_report(self.nifti(AcquisitionTime="10:20:30.500000"), self.dat)
        self.assertTrue(report.startswith("\nWARNING: T1 acquired at 10:20:30, before"))
        self.assertIn("created at 11:00:00", report)

    def test_nifti_acquired_later_or_without_time_is_not_compared(self):
        write_twix(self.dat, LATE)
        for sidecar in (dict(AcquisitionTime="11:20:00.000000"), dict(EchoTime=0.002), {}):
            report = check.frame_report(self.nifti(**sidecar), self.dat)
            self.assertTrue(report.startswith("T1 frame of reference not compared"), report)
            self.assertNotIn("\n", report)

    def test_errors_become_one_line_and_exit_zero(self):
        write_dicom(os.path.join(self.t1, "1.dcm"), "")
        cases = [(self.t1, self.path("missing.dat")), (self.t1, self.path("csi.mat"))]
        write_twix(self.dat, LATE)
        cases += [(self.t1, self.dat), (self.path("empty"), self.dat)]
        os.mkdir(self.path("empty"))
        for t1, dat in cases:
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                self.assertEqual(check.main(["check_t1_frame.py", t1, dat]), 0)
            self.assertTrue(out.getvalue().startswith("T1 frame of reference not checked: "))
            self.assertEqual(out.getvalue().count("\n"), 1)


if __name__ == "__main__":
    unittest.main()
