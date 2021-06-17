function [flatVect] = flatStruct2Vect(flatStruct)
%FLATSTRUCT2VECT 
flatVect = [flatStruct.origin(:); flatStruct.orthProjection(:)];
end

