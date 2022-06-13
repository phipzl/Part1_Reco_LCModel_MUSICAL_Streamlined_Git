function header_rda = convert_header(h)

symbol          = [char(13),char(10)];
header_rda      = cat(2,'>>> Begin of header <<<');

for i = 2:max(size(h))
    header_rda = cat(2,header_rda,symbol,h{i,1},': ',h{i,2});
end

header_rda  = cat(2,header_rda,symbol,'>>> End of header <<<',symbol);