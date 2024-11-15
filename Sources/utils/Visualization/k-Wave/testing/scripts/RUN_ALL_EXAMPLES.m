function RUN_ALL_EXAMPLES
% 11th February 2011
% function to load and run all the examples in the k-Wave Toolbox
% Bradley Treeby

% get examples directory
exam_dir = getkWavePath('examples');

% move there
cd(exam_dir);

% get a list of m-files in the examples directory
filenames = what;
filenames = filenames.m;

% run the examples one by one
for filename_index = 1:length(filenames)

    % save the variables to the workspace to circumvent clear all
    save('RUN_ALL_EXAMPLES_VARS');
    
    % trim the .m extension
    fn = filenames{filename_index};
    fn = fn(1:end-2);
    
    % display the filename
    disp(['Running example ' fn]);

    % run the file
    run(fn);
    
    % load back the variables
    load('RUN_ALL_EXAMPLES_VARS');
    
    % close the figure windows
    close all;
end