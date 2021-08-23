% HWWA_ABBREVIATED_PIPELINE
%
%   This script preprocesses behavioral data created by the hww_gng task.
%   It generates a series of so-called 'intermediate files' -- files that
%   contain small component pieces of a complete data set.
%
%   You can use this script as a template -- however, it's best to save a
%   copy of this file and make changes to that one, leaving the original
%   template intact.
%
%   Unless manually configured, this script will load and save data to a
%   folder called 'data' located alongside the 'script' folder in the hwwa
%   repository folder. To get started, place all raw data (.mat and .edf 
%   files obtained from a given run of the task code) in a folder
%   called 'raw', stored in the 'data' subfolder, and then run this script.
%
%   I.e., by default, your folder structure should look like this:
%
%     hwwa/
%       +hwwa/
%       script/
%       data/
%         raw/
%           file1.mat
%           file1.edf
%           ...
%
%   The full path to the 'data' folder is called the root data directory,
%   and is configurable.
%
%   See also hwwa.config.create, hwwa.dataroot, hwwa.set_dataroot,
%     hwwa.config.load

inputs = hwwa.get_common_make_defaults();

% Creates a file from which most other files are derived; loads raw data.
hwwa.make_unified( inputs );

% Saves .edf files as .mat files
hwwa.make_edfs( inputs );

% Saves task events (in matlab time) as a matrix.
hwwa.make_events( inputs );

% Expresses matlab event times in terms of eyelink's clock.
hwwa.make_alternate_el_events( inputs );

% Converts trial-data to categorical matrix.
hwwa.make_labels( inputs );

% Aligns eye samples to events, creating MxN matrices of M-trials by
% N-samples. In this case, matrices are only generated for
% 'go_target_onset' and 'go_target_acquired' events.
res = hwwa.make_edf_trials( inputs ...
  , 'event',        {'go_target_onset', 'go_target_acquired'} ...
  , 'event_subdir', 'edf_events' ...  % Pull events from 'edf_events' intermediate folder
  , 'look_back',    0 ... % ms
  , 'look_ahead',   1000 ...  % ms
);