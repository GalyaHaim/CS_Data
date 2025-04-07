% Generates a Sensing Matrix for compressive sensing 
% INPUT
% 'prob_size' - how many points in final reconstruction (integer)
% 'num_projs' - how many measurements (projections) you are taking (integer)
% 'type' - describes the elements of 'sensing_mtx' must be one of
% following:
%   (1) 'Hadamard' - Randomly permuted Hadamard matrix (prob_size must be
%          power of 2 unless 'Non-Dense' option is chosen)
%   (2) 'Gaussian' - followed by a number indicating 'sigma' of normal
%        distribution
%   (3) 'Spikes' - Equal Amplitudes, already non-dense
%   (4) 'Poisson'- Spikes that are poisson distributed in amplitude. option
%       must be followed by a number indicating 'lambda' of distribution
%   (5) 'On/Off Spikes' - Randomly chosen even  # indices with half on
%       and half off.
%
% OPTIONAL 
% If a structure generated as output by a previous run is used ALL
% ARGUMENTS are ignored and the same sensing matrix is generated. 
%
% generateSensingMtx(...,'verbose',bool,...) outputs (or suppresses) warnings and other helpful
% information to terminal. Default bool = true
%
% generateSensingMtx(...,'Non-Dense',value,...) generates a sensing matrix with
% guaranteed zeros in each row accoridng to 'value'. 
%       if 'value' is an integer there will be at least that number of
%       nonzeros in each row
%       if 'value' is a float between zero and one at least that fraction
%       of the total number of elements in each row will be zero
%       This option is ignored if 'type' is 'Spikes'.
%       If 'type' is hadamard 'arg' is auto adjusted to ensure that the 
%       number of nonzero elements is a power of 2.
%
% generateSensingMtx(...,'Seed',value,..) uses value as a seed for the 
% random number generator by calling 'rng(value)' for deterministic behavior
% during a single session of MATLAB.
%
% OUTPUT
% 'sensing_mtx' - matrix of size [num_projs,prob_size] with elements chosen
% according to 'type' and 'varargin'
%
% 'params' - structure containing all information about 'sensing_mtx'. If
% you want to generate the same sensing matrix call generateSensingMtx with
% 'params' as the first optional argument. 
function [sensing_mtx, mtx_params] = generateSamplingMtx(prob_size, num_projs, type, varargin)
    mtx_params.rng_settings = rng;
    mtx_params.type = type;
    mtx_params.prob_size = prob_size;
    mtx_params.num_projs = num_projs;
    is_verbose = true;
    nonDenseFlag = false;
    weight_flag = false;
    for i=1:length(varargin)
        if ischar(varargin{i})
            if strcmp(varargin{i},'verbose')
                is_verbose = varargin{i+1};
            end
        end
    end
    %% Check if a structure was passed to varargin and use its parameters
    if ~isempty(varargin)
        if isstruct(varargin{1})
            if is_verbose
                disp('All arguments overwritten, generating Sensing Matrix from a previous result.')
            end
            use_prev = true;
            prob_size = varargin{1}.prob_size;
            num_projs = varargin{1}.num_projs;
            nzeros = varargin{1}.nzeros;
            if nzeros~=0
               nonDenseFlag=true; 
            end
            rng(varargin{1}.rng_settings)
            type = varargin{1}.type;
            switch type
                case 'Gaussian'
                    sigma = varargin{1}.sigma;
                case 'Spikes'
                    nspikes = varargin{1}.nspikes;
                case 'Poissonian'
                    lambda = varargin{1}.lambda;
            end
        else
            %% Check if Sensing Matrix will be Non-Dense and set Random Seeding
            % check if the sensing matrix will have zero elements in projection
            % 'Non-Dense'
            % if Sensing Matrix will be non-dense the next argument should be
            % an integer or a fraction between zero and one indicating what number
            % of elements will be zero in each projection.
            %
            % For deterministic behavior can pass the string 'Seed' followed by
            % your desired settings for MATLAB's 'rng'
            use_prev = false;
            for i=1:length(varargin)
                if ischar(varargin{i})
                    if strcmp(varargin{i},'Non-Dense')
                        if (strcmp(type,'Spikes')==true) || (strcmp(type,'On/Off Spikes')==true)
                            if is_verbose
                                disp('Spike Basis is already Non-Dense. Ignoring `Non-Dense` option.')
                            end
                        else
                            nonDenseFlag = true;
                            density = varargin{i+1};
                            if floor(density)==density
                                nzeros = prob_size-density;
                            elseif isnumeric(density)
                                if abs(density-0.5)<0.5
                                    nzeros = ceil(density*prob_size);
                                else
                                    error('Non-Dense should be followed by an integer or a float between 0 and 1.')
                                end
                            else
                                error('Non-Dense should be followed by an integer or a float between 0 and 1.')
                            end
                        end
                    elseif strcmp(varargin{i},'Seed')
                        rng(varargin{i+1});
                        mtx_params.rng_settings = rng;
                    elseif strcmp(varargin{i}, 'Weights')
                        weights = varargin{i+1};
                        if length(weights)~=prob_size
                            error('Weight vector should be same size as problem size.')             
                        end
                        mtx_params.weights = weights;
                        weight_flag = true;
                    end
                end
            end
            
        end
    end
    %% save Matrix Parameters so far
    if nonDenseFlag
        mtx_params.nzeros = nzeros;
    else
        nzeros = 0;
        mtx_params.nzeros = 0;
    end
    %% Generate Sensing Matrix According to 'type'
    switch type
        case 'Hadamard'
            if nzeros==0
                if floor(log2(prob_size))~=log2(prob_size)
                    error('Hadamard must be a power of 2.')
                end
                %Generate Randomly Permuted Hadamard
                num_nonzero = prob_size;
                % Generate Randomly Permuted Hadamard
                %choose random rows ignoring the first row
                Q=randperm(num_nonzero);
                had=seq2had(num_nonzero);
                m=Q(1:num_projs); %rows to take
                for i=1:num_projs
                    if Q(i)==1
                        m(i)=Q(num_projs+1); %ignore first row if it occurs
                    end
                end
                clear Q
                p=randperm(num_nonzero); %random permutation matrix for switching rows
                sensing_mtx = zeros(num_projs,  num_nonzero);
                sensing_mtx(:,p)=mRowHn(had(m), num_nonzero);
                
            else
                
                % Change number of zero elements per projection to ensure
                % Hadamard is a power of 2.
                if is_verbose
                    disp('WARNING: Changing number of zero elements so problem size is power of 2 for Hadamard.')
                end
                
                if weight_flag
                    if is_verbose
                        disp('Weights given.')
                    end
                    num_nonzero = 2^(floor(log2(prob_size-nzeros)));
                    mtx_params.nzeros = prob_size-num_nonzero;
                    had=seq2had(num_nonzero);                  
                    sensing_mtx = zeros(num_projs,  prob_size);
                    for i=1:num_projs
                        m=randi(num_nonzero);
                        while m==1
                            m=randi(num_nonzero);
                        end
                        weighted_indices = datasample(1:prob_size,num_nonzero,'Replace',false,'Weights',weights);
                        sensing_mtx(i,weighted_indices)=mRowHn(had(m), num_nonzero);
                    end                
                else
                    num_nonzero = 2^(floor(log2(prob_size-nzeros)));
                    mtx_params.nzeros = prob_size-num_nonzero;
                    
                    had=seq2had(num_nonzero);
                    sensing_mtx = zeros(num_projs,  prob_size);
                    for i=1:num_projs
                        m=randi(num_nonzero);
                        while m==1
                            m=randi(num_nonzero);
                        end
                        p=randperm(prob_size); %random permutation matrix for switching rows
                        sensing_mtx(i,p(1:num_nonzero))=mRowHn(had(m), num_nonzero);
                    end
                end
                
            end
            
        case 'Gaussian'
            if ~use_prev
                sigma=varargin{1};
                mtx_params.sigma = sigma;
            end
            if nzeros==0
                sensing_mtx = sigma*randn(num_projs, prob_size);
            else
                sensing_mtx=zeros(num_projs,prob_size);
                num_nonzero = prob_size-nzeros; %number of nonzero elements in each row
                for i=1:num_projs
                    picks = randperm(prob_size);
                    sensing_mtx(i,picks(1:num_nonzero))=sigma*randn(1,num_nonzero);
                end
            end
            
        case 'Spikes'
            if ~use_prev
                nspikes=varargin{1};
                mtx_params.nspikes = nspikes;
            end
            sensing_mtx=zeros(num_projs,prob_size);
            for i=1:num_projs
                picks = randperm(prob_size);
                sensing_mtx(i,picks(1:nspikes))=1;
            end
            
        case 'Poissonian'
            % one parameter 'lambda' is needed for poisson distribution
            if ~use_prev
                lambda=varargin{1};
                mtx_params.lambda = lambda;
            end
           
            if nzeros==0
                num_nonzero = prob_size;
                sensing_mtx = poisspdf(poissrnd(lambda,num_projs,num_nonzero),lambda);
            else
                sensing_mtx=zeros(num_projs,prob_size);
                num_nonzero = prob_size-nzeros; %number of nonzero elements in each row
                for i=1:num_projs
                    picks = randperm(prob_size);
                    sensing_mtx(i,picks(1:num_nonzero))=poisspdf(poissrnd(lambda,1,num_nonzero),lambda);
                end
            end
        case 'On/Off Spikes'
            % On/Off Spikes
            if ~use_prev
                nspikes=varargin{1};
                mtx_params.nspikes = nspikes;
                num_nonzero = nspikes;
            end
            mtx_params.nzeros = prob_size-num_nonzero;
            sensing_mtx = zeros(num_projs,  prob_size);
            num_pos_indices = num_nonzero;
            if weight_flag
                if is_verbose
                    disp('Weights given.')
                end
                for i=1:num_projs
                    weighted_indices = datasample(1:prob_size,num_nonzero,'Replace',false,'Weights',weights);
                    pos_indices = datasample(weighted_indices,num_pos_indices,'Replace',false);
                    sensing_mtx(i,pos_indices) = 1;
                    sensing_mtx(i,setdiff(weighted_indices,pos_indices)) = 0;
                end                
            else
                for i=1:num_projs
                    weighted_indices = datasample(1:prob_size,num_nonzero,'Replace',false);
                    pos_indices = datasample(weighted_indices,num_pos_indices,'Replace',false);
                    sensing_mtx(i,pos_indices) = 1;
                    sensing_mtx(i,setdiff(weighted_indices,pos_indices)) = 0;
                end                            
            end
           
        otherwise
            error(['Unknown Sensing Matrix Type, valid choices: `Hadamard`',...
                ', `Gaussian`, `Spike`, `Poissonian`, `On/Off Spikes`'])
    end