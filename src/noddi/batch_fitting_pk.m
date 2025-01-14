function batch_fitting_pk(roifile, bval_vol, bvec_vol, model, outputfile, poolsize)
%
% function batch_fitting(roifile, protocol, model, outputfile, poolsize)
%
% This function does batch fitting to the voxels in an entire ROI created
% with CreateROI function.
%
% This function supports parallel computation on multiple CPUs by taking
% advantage of Matlab's Parallel Computing Toolbox; the availability of
% this toolbox is thus a requirement. If you don't have access to this
% toolbox, please use batch_fitting_single instead.
%
% The function can handle interrupted runs.  If a run is stopped
% prematurely, either through an unexpected sytem crash or through a
% deliberate interruption, simply re-run the same command to resume the
% fitting.
%
%
% Input:
%
% roifile: the ROI file created with CreateROI
%
% protocol: the protocol object created with FSL2Protocol
%
% model: the model object created with MakeModel
%
% outputfile: the name of the mat file to store the fitted parameters
%
% poolsize (optional): the number of parallel processes to run
%
%
% author: Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%

% first check if there is a file there to resume
if exist(outputfile, 'file')
    output = load(outputfile);
    if isfield(output, 'split_end')
        % previously run and need to be restarted
        split_end = output.split_end;
        gsps = output.gsps;
        mlps = output.mlps;
        fobj_gs = output.fobj_gs;
        fobj_ml = output.fobj_ml;
        error_code = output.error_code;
        if model.noOfStages == 3
            mcmcps = output.mcmcps;
        end
        current_split_start = split_end + 1;
        fprintf('Resume an interrupted run from %i\n', current_split_start);
    else
        % completed
        fprintf('An output file of the same name detected.\n');
        fprintf('Choose a different output file name.\n');
        return;
    end
else
    % if this is the first run
    current_split_start = 1;
end

% initiate the parallel environment if necessary
pool = gcp('nocreate');
if isempty(pool)
    if (nargin < 5)
        pool = parpool('local');
    else
        pool = parpool('local', poolsize);
    end
end

% load the roi file
input = load(roifile);
roi = input.roi;
indx = input.idx;
numOfVoxels = size(roi,1);

% set up the fitting parameter variables if it is the first run
if current_split_start == 1
    gsps = zeros(numOfVoxels, model.numParams);
    mlps = zeros(numOfVoxels, model.numParams);
    fobj_gs = zeros(numOfVoxels, 1);
    fobj_ml = zeros(numOfVoxels, 1);
    error_code = zeros(numOfVoxels, 1);
    if model.noOfStages == 3
        mcmcps = zeros(numOfVoxels, model.MCMC.samples, model.numParams + 1);
    end
end

% set up the PARFOR Progress Monitor
progressStepSize = 100;
ppm = ParforProgMon(['Fitting ' roifile, ' : '], numOfVoxels-current_split_start+1,...
                    progressStepSize, 400, 80);

tic

fprintf('%i of voxels to fit\n', numOfVoxels-current_split_start+1);

% start the parallel fitting
for split_start=current_split_start:progressStepSize:numOfVoxels
    % set up the split end
    split_end = split_start + progressStepSize - 1;
    if split_end > numOfVoxels
        split_end = numOfVoxels;
    end
    
    % fit the split
    parfor i=split_start:split_end
        
        % get the MR signals for the voxel i
        voxel = roi(i,:)';
        voxel_idx_x = indx(i,1);
        voxel_idx_y = indx(i,2);
        voxel_idx_z = indx(i,3);
        
        % get the protocol for that voxel
        bval_voxel = bval_vol(voxel_idx_x,voxel_idx_y,voxel_idx_z,:);
        bvec_voxel = bvec_vol(voxel_idx_x,voxel_idx_y,voxel_idx_z,:,:);
        protocol = FSL2Protocol_pk(bval_voxel, bvec_voxel);
        
        % fit the voxel
        if model.noOfStages == 2
            [gsps(i,:), fobj_gs(i), mlps(i,:), fobj_ml(i), error_code(i)] = ThreeStageFittingVoxel(voxel, protocol, model);
        else
            [gsps(i,:), fobj_gs(i), mlps(i,:), fobj_ml(i), error_code(i), mcmcps(i,:,:)] = ThreeStageFittingVoxel(voxel, protocol, model);
        end
        
        % report to the progress monitor
        if mod(i, progressStepSize)==0
            ppm.increment();
        end
        
    end
    
    % save the temporary results of the split
    if model.noOfStages == 2
        save(outputfile, 'split_end', 'model', 'gsps', 'fobj_gs', 'mlps', 'fobj_ml', 'error_code');
    else
        save(outputfile, 'split_end', 'model', 'gsps', 'fobj_gs', 'mlps', 'fobj_ml', 'mcmcps', 'error_code');
    end
    
end

toc

ppm.delete();

% save the fitted parameters
if model.noOfStages == 2
    save(outputfile, 'model', 'gsps', 'fobj_gs', 'mlps', 'fobj_ml', 'error_code');
else
    save(outputfile, 'model', 'gsps', 'fobj_gs', 'mlps', 'fobj_ml', 'mcmcps', 'error_code');
end

% close the parallel pool
delete(pool);

