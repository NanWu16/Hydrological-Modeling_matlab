function varargout = ini_zero_2dim(Nx, Ny, Number)
    for i = 1 : Number
        eval(strcat("a_", num2str(i), " = zeros(Nx, Ny);"))
        eval(strcat("varargout{i} = a_", num2str(i), ";"))
    end
end