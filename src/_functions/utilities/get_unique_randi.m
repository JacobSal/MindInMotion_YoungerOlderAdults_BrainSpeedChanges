function [a1] = get_unique_randi(varargin)
%GET_UNIQUE_RANDI wrapper for randi.m matlab function that generates unique
%integers for a given randi input. 
% varargin = ...
% (imax)
% (imax,n)
% (imax,n,sz1,...,szN)
% (imax,sz)
% (imax,...,classname)
tic

a1 = randi(varargin{:});
b1 = unique(a1);
while length(a1) ~= length(b1)
    a1 = randi(varargin{:});
    b1 = unique(a1);
end
toc
end

