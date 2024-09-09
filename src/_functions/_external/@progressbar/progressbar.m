function p=progressbar(pr)
if nargin==0
    p.progress_bar_position=0;
    p.relapsed_time=0.01;
    p=class(p,'progressbar');
    tic;
elseif isa(pr,'progressbar')
   p=pr;
end
end