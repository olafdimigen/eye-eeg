% PARSEQ - parses a column vector according to its incremental sequence of 
% unchanging values. Use for identifying sequential blocks of data. 
%
% Inputs: 
%   a:  column vector
% Outputs:
%   b1: begin 
%   b2: end 
%   b3: length of data block
%
% Reinhold Kliegl & Jochen Laubrock, 04-10-2004

function b = parseq(a)
tmp=find(diff(a)~=0);
b=[[1; tmp+1] [tmp; length(a)]];
b=[b b(:,2)-b(:,1)+1];