% findsequence2(): finds continuous sequences of integer values in
% monotonically increasing indices.
% Returns begin, end and length of sequences
%
% Author: ur, reinacul@psychologie.hu-berlin.de
%
% Input:
%   x:        column vector of integer values, e.g. x = [1;2;3;4;7;8]
%
% Outputs:
%   out(:,1): first value of continuous sequence (not: first index)
%   out(:,2): last value of continuous sequence (not: last index)
%   out(:,3): length of continuous sequence

function [out] = findsequence2(x)

if size(x,2)>1
    if size(x,1)>1
        error('findsequence2() expects a vector as input')
    else
        fprintf('\nFunction findsequence2() expects column vector as input');
        fprintf('\nInput is transposed');
        x = x';
    end
end

indBegin = [1;1+find(diff(x)~=1)];
indEnd   = [find(diff(x)~=1);length(x)];
lenSeq   = (indEnd-indBegin)+1;
out      = [x(indBegin),x(indEnd),lenSeq];