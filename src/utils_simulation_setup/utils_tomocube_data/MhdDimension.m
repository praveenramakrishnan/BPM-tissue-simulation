function [out] = MhdDimension(mhd)
    fp = fopen(mhd);
    tline = fgetl(fp);
    while(ischar(tline))
        [tok, remain] = strtok(tline, ' = ');
        if length(tok)==7
            if tok == 'DimSize'
                remain = remain(4:end);
                idx = 1;
                while ~isempty(remain)
                    [tok, remain] = strtok(remain, ' '); 
                    out(idx) = str2num(tok);
                    idx = idx + 1;
                end
            end
        end
        tline = fgetl(fp);
    end
    fclose(fp);
end