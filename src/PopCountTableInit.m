function PopCountStruct = PopCountTableInit(TableShift)
PopCountStruct = struct('TShift',TableShift,'TSize',[],'Table',[]);
PopCountStruct.TSize = bitshift(uint32(1),TableShift);
PopCountStruct.Table = zeros(PopCountStruct.TSize,1);
for i = 0:PopCountStruct.TSize
    PopCountStruct.Table(i+1) = PopCountReal(i);
end
end

function counts = PopCountReal(w)
binary = dec2bin(w);
counts = count(binary,'1');
end