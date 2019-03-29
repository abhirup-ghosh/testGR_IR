% Matlab script to generate documentation with m2html

addpath ./m2html/
m2html('mfiles','EOB', 'htmldir','Doc', 'recursive','on', 'global','on');
