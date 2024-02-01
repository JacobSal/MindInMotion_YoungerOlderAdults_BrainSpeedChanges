function [out] = chk_cell(x,cl_i)
    if x == cl_i
        out = 1;
    else
        out = 0;
    end
    out = logical(out);
end