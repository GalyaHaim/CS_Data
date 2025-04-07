classdef spike_adaptive_l1homotopy
    %ADAPTIVE_RECONSTRUCTION Summary of this class goes here
    %   Detailed explanation goes here
    %latest-used with simulations
    properties
        smp_mtx {mustBeNumeric}
        smp_freqs (:,1) {mustBeNumeric}
        
        
        
        % BASIS
        dual_spike_basis = false;
        mdpt;
        basis_mtx {mustBeNumeric}
        basis_freqs (:,1) {mustBeNumeric}
        
        sparse_rep (:,1) {mustBeNumeric}
        
        projs (:,1) {mustBeNumeric}
        measurement_number {mustBeInteger} = 0
        Continue {mustBeNumericOrLogical} = true
        est_mean (1,1) {mustBeNumeric} = nan
        est_std (1,1) {mustBeNumeric} = nan
        mean_confidence (1,1) {mustBeNumeric} = Inf;
        
        reconstruction (:,1) {mustBeNumeric}
        num_sample_freqs (1,1) {mustBeInteger} = 1;
        guess_contrast = false;
        contrast (1,1) {mustBeNumeric} = 0.0008; %30/6
        found_fit = false;
        
        % final fit data
        fit (:,1) {mustBeNumeric}
        pk_locs (:,1) {mustBeNumeric}
        pk_unc (:,1) {mustBeNumeric}
        width (:,1) {mustBeNumeric}
        width_unc (:,1) {mustBeNumeric}
        
        ref_measurements (:,1) {mustBeNumeric}
        sig_measurements (:,1) {mustBeNumeric}
        
        min_height_window (1,1) {mustBeNumeric} = 0.00;
        all_pks = [];
        all_pks_unc = [];
        % did the 4 previous fits all find the same peaks
        fit_converged = false;
        % how many times did we successfully fit every step
        all_pass = 0;
        
        %target error (MHz)
        pk_target_error = 5;
        
        % if the correct number of peaks is not found in a step
        % we can save them anyway.
        candidate_locations = [];
        first_pass_peaks = [];
        second_pass_peaks = [];
        
        
        % for last step when we fit widths
        basis_widths (:,1) {mustBeNumeric}
        
        % for saving every step
        all_reconstructions = [];
        
        % when to start saving statistics
        start_statistics = nan;
        save_statistics = false;
        
        
        %% FOR SIMULATION DATA
        sim
        get_err = false;
        actual_err = [];
        actual_err_locations = [];
        meas_num = [];
        all_err = [];

    end
    
    properties (Hidden)
        % whether to plot at each step
        display = false;
        % stores random state in case you want to repeat experiment with
        % deterministic behavior
        rng_state = rng;
        % default store entire sensing matrix and basis in memory
        % otherwise will calculate on fly
        large_scale = false;
        % sampling weights
        weights (:,1) {mustBeNumeric}
        % ratio to start first reconstruction
        min_ratio = 0.07;
        % swithces to true after min_ratio is reached.
        flag_start_reconstructions = false;
        % ratio to end reconstruction
        max_ratio = 0.25;
        % how many frequencies we can measure at
        num_pts = 0;
        % number of intial projections before we start reconstructions
        num_initial_projs (1,1) {mustBeNumeric}
        % number of intial projections before we start reconstructions
        num_max_projs (1,1) {mustBeNumeric}
        % lorentzian kernel
        basis_kernel
        % number of lorentzians in initial guess
        problem_size = 651;
        % width initial guess
        width_initial = 12;
        % array of guessed windows
        windows
        % number of good guesses for a window occurs
        good_guess_flag (1,1) {mustBeInteger} = 0;
        % number of zeros that appear in each projection
        num_zeros = 0;
        % flag whether window is good enough to start looking for peaks
        look_for_peaks = false;
        target_num_pks = 8
        %window around each candidate peak
        freq_window = 10;
        num_bins = 20;
        % for fitting widths around each peak
        %         width_window = 5;
        %         num_width_bins = 50;
        width_window = 2;
        num_width_bins = 20;
        
        % keep track of means & variances
        all_means (:,1)
        all_stds (:,1)
        all_contrasts (:,1)
        all_confidences (:,1)
        
        peaks_found = false;
        final_window
        
        % for plotting
        
        fig
        ax
        
        
        %whether to use candidate peak locations to tighten window
        
        use_peaks = true;
        poll_peaks_randomly = false;
        poll_weakest_peak = false;
        poll_window = false;
        stop_polling_peaks = false;
        polled_peaks = [];
        curr_locs = [];
        curr_locs_unc = [];
        curr_widths = [];
        curr_widths_unc = [];
        curr_fit = [];
        found_candidate_peaks = false;
        fit_widths = false;
        num_second_pass = 0;
        
        
        % if we want to output a movie
        get_movie = false;
        animation = struct('cdata',[],'colormap',[]);
        video
        
        % midpoints
        midpoints = [];
        
        window_start
        window_end
        
        %for NESTA
        opts
        AA
        AAt
        
        % For L1 Homotopy
        in
        
        indice = [];
        
        sig_mean = 0;
        ref_mean = 0;
        
        % freq window
        
        
        % anonymous function from 'lorentzian fit'
        L_func
    end
    
    methods
        function exp = SaveExperiment(obj,fileName)
            % Save all properties as exp.propertyName, and adds a exp.date
            % field. If file name is given - sames using filename. Elese -
            % open a GUI interface
            % saves using
            names = fields(obj);
            for k = 1:length(names)
                command = sprintf('exp.%s = obj.%s;',names{k},names{k});
                eval(command);
            end
            exp.date = datestr(now,'yyyy-mm-dd_HH-MM-SS');
            eval(sprintf('save(''%s'',''exp'');',fileName));
        end
        function obj = spike_adaptive_l1homotopy(freqs, num_sample_freqs, contrast, window_start, window_end, varargin)
            %ADAPTIVE_RECONSTRUCTION Construct an instance of this class
            %   If no argument is given runs a simulation based on
            %   previously taken data. Otherwise it takes exactly one
            %   argument of the vector of frequencies we will be
            %   measuring at
            
            % get frequencies, only requirement is that they are
            % unique and all the peaks should lay between the minimum
            % and maximum
            if any(diff(freqs)==0)
                error('Frequencies should be unique.')
            end
            obj.smp_freqs = freqs; %freq_range
            obj.num_pts = length(freqs);
            obj.num_sample_freqs = num_sample_freqs;
            obj.contrast = contrast;
            obj.window_start = window_start;
            obj.window_end = window_end;
            
            % get initial sampling matrix
            obj.num_initial_projs = ceil(obj.min_ratio*obj.num_pts);
            obj.num_zeros = obj.num_pts-obj.num_sample_freqs;
            
            [obj.smp_mtx, ~] = generateSamplingMtx(obj.num_pts, obj.num_initial_projs, 'On/Off Spikes',obj.num_sample_freqs,...
                'verbose',false);
            
            % initial basis matrix
            obj.basis_kernel =  lorentzian_kernel('ZMN');
            
            obj.basis_freqs = linspace(window_start, window_end, obj.problem_size);
            
            obj.basis_mtx = get_basis_mtx(obj.basis_kernel, obj.smp_freqs, obj.basis_freqs, ...
                obj.width_initial);
            
            % failsafe stopping condition
            
            obj.num_max_projs = ceil(obj.max_ratio*obj.num_pts);
            if obj.num_max_projs<100
                if 0.5*obj.num_pts>100
                    obj.num_max_projs = 100;
                else
                    obj.num_max_projs = ceil(0.5*obj.num_pts);
                end
                
            end
            
            
            % TVAL3 opts
            obj.opts.mu = 2^9;
            obj.opts.beta = 2^6;
            obj.opts.mu0 = 2^3;       % trigger continuation shceme
            obj.opts.beta0 = 2^-2;    % trigger continuation shceme
            obj.opts.maxcnt = 300;
            obj.opts.tol_inn = 1e-3;
            obj.opts.tol = 1E-5;
            obj.opts.maxit = 1000;
            obj.opts.TVnorm = 1;
            obj.opts.TVL2 = true;
            obj.opts.nonneg = true;
            obj.opts.scale_b = false;
            obj.opts.scale_A = false;
            obj.opts.disp = false;
            
            % L1 Homotopy Option
            obj.in = [];
            
            % final fit data
            obj.fit = zeros(obj.num_pts,1);
            obj.pk_locs = zeros(obj.target_num_pks,1);
            obj.pk_unc = zeros(obj.target_num_pks,1);
            obj.width = zeros(obj.target_num_pks,1);
            obj.width_unc = zeros(obj.target_num_pks,1);
            obj.all_means = zeros(obj.num_max_projs,1);
            
            % figure for plotting too
            if obj.display
                obj.fig = figure;
                if obj.get_movie
                    axis tight manual
                    xlim([obj.smp_freqs(1) obj.smp_freqs(end)])
                    ylim([-0.02 0.055])
                    obj.ax = gca;
                    obj.ax.NextPlot = 'replaceChildren';
                    obj.video = VideoWriter('Compressed_Sensing','MPEG-4');
                    obj.video.FrameRate = 5;
                    obj.video.Quality = 100;
                    %                 obj.video.LosslessCompression = true;
                    open(obj.video)
                else
                    obj.ax = gca;
                end
                
            end
            
            obj.all_pks = [];
            obj.all_pks_unc = [];
            
            
            %% If there is a passed argument assume it is a simulation
            if ~isempty(varargin)
                obj.sim = varargin{1};
                obj.get_err = false;
            end
        end
        %%
        function [obj, freqs_to_measure] = frequency(obj)
            % send next frequencies to measure at
            freqs_to_measure = zeros(obj.num_sample_freqs,1);
            if (~obj.flag_start_reconstructions) && obj.Continue
                % we haven't measured enough points to start reconstructing
                % yet.
                obj.measurement_number=obj.measurement_number+1;
                freqs_to_measure(1:obj.num_sample_freqs) = ...
                    obj.smp_freqs(obj.smp_mtx(obj.measurement_number,:)>0);
                if obj.measurement_number==obj.num_initial_projs
                    obj.flag_start_reconstructions=true;
                end
            elseif obj.Continue
                if obj.use_peaks && obj.poll_peaks_randomly
                    %% TODO: use peak locations
                    % need to make this only do one peak
                    num_pks = length(obj.curr_locs);
                    obj.windows = getWindows(obj.smp_freqs,...
                        obj.curr_locs(randi(num_pks)),10);
                    if ~isempty(obj.polled_peaks)
                        while sum(obj.polled_peaks.*obj.windows)~=0
                            obj.windows = getWindows(obj.smp_freqs,...
                                obj.curr_locs(randi(num_pks)),10);
                        end
                    end
                    obj.polled_peaks = obj.windows;
                    
                    windowed_row =  generateSamplingMtx(obj.num_pts, 1, ...
                        'On/Off Spikes',1,'verbose',false,...
                        'Weights',obj.windows);
                    
                    
                    if obj.num_sample_freqs==1
                        new_row = windowed_row;
                    else
                        obj.weights = ones(obj.num_pts,1);
                        obj.weights = double(xor(obj.weights, obj.windows));
                        complement_row =   generateSamplingMtx(obj.num_pts, 1, ...
                            'On/Off Spikes',obj.num_sample_freqs-1,'verbose',false,...
                            'Weights',obj.weights);
                        new_row = complement_row+windowed_row;
                    end
                    obj.smp_mtx = [obj.smp_mtx; new_row];  % add to sampling matrix
                    obj.measurement_number=obj.measurement_number+1;
                    freqs_to_measure(1:obj.num_sample_freqs) = ...
                        obj.smp_freqs(obj.smp_mtx(obj.measurement_number,:)>0);
                    
                    
                elseif obj.poll_window
                    midpoint = 0.5*(obj.window_end+obj.window_start);
                    window_width = obj.window_end-obj.window_start;
                    obj.windows = getWindows(obj.smp_freqs,midpoint,window_width);
                    occurs = 1./(1+(getOccurs(obj.smp_mtx)').^4);
                    obj.windows = obj.windows.*occurs;
                    new_row =  generateSamplingMtx(obj.num_pts, 1, ...
                        'On/Off Spikes',1,'verbose',false,...
                        'Weights',obj.windows);
                    
                    if obj.num_sample_freqs>1
                        obj.windows = getWindows(obj.smp_freqs,obj.smp_freqs(logical(new_row)),obj.freq_window);
                        obj.weights = ones(obj.num_pts,1);
                        obj.weights = double(xor(obj.weights, obj.windows));
                        complement_row =   generateSamplingMtx(obj.num_pts, 1, ...
                            'On/Off Spikes',obj.num_sample_freqs-1,'verbose',false,...
                            'Weights',obj.weights);
                        new_row = complement_row+new_row;
                    end
                    obj.smp_mtx = [obj.smp_mtx; new_row];  % add to sampling matrix
                    obj.measurement_number=obj.measurement_number+1;
                    freqs_to_measure(1:obj.num_sample_freqs) = ...
                        obj.smp_freqs(obj.smp_mtx(obj.measurement_number,:)>0);
                else
                    obj.weights = ones(obj.num_pts,1);
                    new_row =  generateSamplingMtx(obj.num_pts, 1, ...
                        'On/Off Spikes',obj.num_sample_freqs,'verbose',false,...
                        'Weights',obj.weights);
                    obj.smp_mtx = [obj.smp_mtx; new_row];  % add to sampling matrix
                    obj.measurement_number=obj.measurement_number+1;
                    freqs_to_measure(1:obj.num_sample_freqs) = ...
                        obj.smp_freqs(obj.smp_mtx(obj.measurement_number,:)>0);
                end
            end
            
            if obj.measurement_number>=obj.num_max_projs
                obj.Continue=false;
            end
        end
        %%
        function obj = process_data(obj, measurement_data)
            if obj.measurement_number>=obj.num_max_projs
                obj.Continue=false;
            end
            % measurement data is passed as an array (y,x,4,1)
            % (y,x) is row and column of a pixel
            % the third dimension is the [positive signal, ...
            % positive_reference] at
            % the point (y,x).
            % currently set-up so that we look over the whole frame
%             [nrows, ncols, ~] = size(measurement_data);
%             pos_proj_sig = squeeze(sum(sum(measurement_data(:,:,1),...
%                 1),2))/(nrows*ncols);
%             pos_proj_ref = squeeze(sum(sum(measurement_data(:,:,2),...
%                 1),2))/(nrows*ncols);

            pos_proj_sig=measurement_data(1);
            pos_proj_ref=measurement_data(2);
            %obj.ref_measurements(obj.measurement_number) =  pos_proj_ref/obj.num_sample_freqs;
            obj.ref_measurements(obj.measurement_number) =  pos_proj_ref;
            obj.ref_mean = mean(obj.ref_measurements);
            
            %obj.sig_measurements(obj.measurement_number) =  pos_proj_sig/obj.num_sample_freqs;
            obj.sig_measurements(obj.measurement_number) =  pos_proj_sig;
            obj.sig_mean = mean(obj.sig_measurements);
            
            % THIS LINE
            obj.projs = obj.projs*obj.est_mean;
            
            
            %obj.windows = zeros(size(obj.smp_freqs));
            % new mean
            %obj.est_mean = (old_mean + pos_proj_ref)./(obj.measurement_number*obj.num_sample_freqs);
            obj.est_mean = mean(obj.ref_measurements);
            
            % new std deviation
            obj.est_std = std(obj.ref_measurements);
            
            % confidence interval in mean
            obj.mean_confidence = 2*tinv(0.99,obj.measurement_number-1)*...
                obj.est_std/(sqrt(obj.measurement_number));
            %obj.mean_confidence = tinv(0.99,obj.measurement_number-1);
            
            
            if obj.guess_contrast
                obj.contrast = obj.est_std./obj.est_mean;
%                 obj.contrast= 
            end
            
            % new projections
            obj.projs(obj.measurement_number) = (obj.num_sample_freqs*(pos_proj_ref-pos_proj_sig));
            %obj.projs(obj.measurement_number) = 1-(pos_proj_sig/pos_proj_ref);
            % THIS LINE
            obj.projs = obj.projs/obj.est_mean;
            
            obj.all_means(obj.measurement_number) = obj.est_mean;
            obj.all_stds(obj.measurement_number) = obj.est_std;
            obj.all_contrasts(obj.measurement_number) = obj.contrast;
            obj.all_confidences(obj.measurement_number) = obj.mean_confidence;
            
            if (obj.flag_start_reconstructions) && (~obj.fit_converged || mod(obj.num_second_pass,3)==0)
                % reconstruct
                %obj.opts.init = 1;
                
                if false && obj.found_candidate_peaks && mod(obj.num_second_pass,3)~=0
                    num_pks = length(obj.curr_locs);
                    detunings = zeros(obj.num_bins*num_pks,1);
                    obj.basis_mtx = zeros(obj.num_pts, obj.num_bins*num_pks);
                    for width_idx=1:num_pks
                        curr_indices = ((width_idx-1)*obj.num_bins+1):(obj.num_bins*width_idx);
                        detunings(curr_indices)=linspace(obj.curr_locs(width_idx)-0.5*obj.freq_window,obj.curr_locs(width_idx)+0.5*obj.freq_window,obj.num_bins);
                        obj.basis_mtx(:,curr_indices) = get_basis_mtx(obj.basis_kernel, obj.smp_freqs, detunings(curr_indices), ...
                            obj.curr_widths(width_idx));
                    end
                else
                    if obj.dual_spike_basis
                        obj.mdpt = 0.5*(obj.basis_freqs(1)+obj.basis_freqs(end));
                        split_freqs = linspace(obj.mdpt,obj.basis_freqs(end),obj.problem_size);
                        obj.basis_mtx = get_dual_basis_mtx(obj.basis_kernel, obj.smp_freqs, split_freqs, ...
                            obj.width_initial, obj.mdpt);
                    else                       
                        obj.basis_mtx = get_basis_mtx(obj.basis_kernel, obj.smp_freqs, obj.basis_freqs, ...
                            obj.width_initial);
                    end
                    
                    
                end
                sensing_mtx = obj.smp_mtx*obj.basis_mtx;
                
                
                
                
                
                %                 delta = 1e-5;
                %                 [U,S,V] = svd(sensing_mtx,'econ');
                %                 obj.opts.USV.U=U;
                %                 obj.opts.USV.S=S;
                %                 obj.opts.USV.V=V;
                %                 obj.AA = @(x) sensing_mtx*x;
                %                 obj.AAt = @(x) sensing_mtx'*x;
                try
                    %                     [obj.sparse_rep, ~] = TVAL3(sensing_mtx,obj.projs-mean(obj.projs),...
                    %                         size(obj.basis_mtx,2),1,obj.opts);
                    %muf = 1e-3;
                    %[obj.sparse_rep,niter,resid,outData] = NESTA(obj.AA,obj.AAt,obj.projs-mean(obj.projs),muf,delta,obj.opts);
                    
                    
                    
                    muf = 1e-3;
                    %tau = 1e-2*max(1e-4*max(abs(sensing_mtx'*(obj.projs))),obj.est_std*sqrt(log(length(obj.smp_freqs*obj.num_sample_freqs))));
                    if obj.get_err && obj.sim.sigma==0
                        sigma = 0.0000001;%00001
%                         sigma = 0.000001;%00001
                    else
                        SNR = (0.01/obj.est_std)/obj.num_sample_freqs;
                        M = obj.measurement_number;
                        sigma = 0.5*sqrt(norm(0.01/8)^2/10^(SNR/10)/M);
                    end                   
                    tau = sigma*sqrt(log(length(obj.smp_freqs)));
                    
                    
                    
                    
                    delx_mode = 'mil';
                    %err_fun = @(z) (norm(x-z)/norm(x))^2;
                    obj.in.tau = tau;
                    obj.in.Te = 8;
                    obj.in.nonneg = 1; % Positivity constraint flag
                    obj.in.delx_mode = delx_mode;
                    obj.in.debias = 1;
                    obj.in.verbose = 0;
                    obj.in.plots = 0;
                    obj.in.record = 0;
                    %in.err_fun = err_fun;
                    x = l1homotopy(sensing_mtx,obj.projs-mean(obj.projs),obj.in);
                    obj.sparse_rep = x.x_out;
                    
                    %
%                     if obj.found_candidate_peaks
%                         obj.opts.init = obj.sparse_rep;
%                     else
%                         obj.opts.init = 1;
%                     end
%                     [obj.sparse_rep, ~] = TVAL3(sensing_mtx,obj.projs-mean(obj.projs),...
%                         size(obj.basis_mtx,2),1,obj.opts);
                    
                    
                    %
                    %                     xh_old = obj.sparse_rep; tau_old = tau;
                    %
                    %
                    %                     W = tau;
                    %                     W_old = tau_old;
                    %                     Atr = sensing_mtx'*(sensing_mtx*obj.sparse_rep-obj.projs-mean(obj.projs));
                    %                     u =  -W.*sign(xh_old)-Atr;
                    %                     pk_old = Atr+u;
                    %
                    %                     obj.in = x;
                    %                     gamma_old = obj.in.gamma;
                    %                     AtAgx = sensing_mtx(:,gamma_old)'*sensing_mtx(:,gamma_old);
                    %                     iAtAgx = inv(AtAgx);
                    %                     obj.in.iAtA = iAtAgx;
                    %
                    %                     obj.in.xh_old = xh_old;
                    %                     obj.in.pk_old = pk_old;
                    %                     obj.in.u = u;
                    %                     obj.in.W = W;
                    %                     obj.in.W_old = W_old;
                    %                     obj.in.delx_mode = delx_mode;
                    %
                    
                catch
                    disp('EW')
                end
                if ~isempty(obj.sparse_rep)
                    obj.reconstruction =  obj.basis_mtx*obj.sparse_rep;
                    
                    if obj.display
                        if sum(obj.fit)==0
                            plot(obj.ax, obj.smp_freqs, obj.reconstruction)
                        else
                            plot(obj.ax, obj.smp_freqs, obj.fit-mean(obj.fit))
                        end
                        title(num2str(obj.measurement_number))
                        drawnow
                    end
                    
                    if obj.measurement_number>=obj.num_max_projs
                        obj.Continue=false;
                        obj.fit = obj.curr_fit;
                        obj.pk_locs = obj.curr_locs;
                        obj.pk_unc = obj.curr_locs_unc;
                        obj.width = obj.curr_widths;
                        obj.width_unc = obj.curr_widths_unc;
                    end
                    % NEED TO PUT IN SOMETHING ABOUT WHETHER FIT LOOKS OK
                    [~, curr_num_pks, curr_guess] = getFitGuess(obj.smp_freqs,... % from the reconstruction takes the initial parameters for the fit
                        obj.reconstruction, obj.contrast);
                    if curr_num_pks == obj.target_num_pks
                        try
                            [obj.curr_fit, curr_params, ~, ~, curr_conf, obj.L_func] = ... %tries fitting on reall data
                                lorentzian_fit_mtx(obj.smp_freqs',...
                                (obj.projs-mean(obj.projs))', 2, 2, curr_num_pks, obj.smp_mtx, ...
                                curr_guess);
                            [obj.curr_locs, obj.curr_locs_unc] = getFitVals(curr_params, ...
                                curr_conf, 'Peak');
                            [obj.curr_widths, obj.curr_widths_unc] = getFitVals(curr_params, ...
                                curr_conf, 'Width');
                            
                            % look at uncertainties and consistency in peaks
                            obj.first_pass_peaks = [obj.first_pass_peaks ; obj.curr_locs];
                            if size(obj.first_pass_peaks,1)>4
                                % look at peak locations over the past 4 fits
                                past_pk_errs = abs(obj.first_pass_peaks(end,:)...
                                    -mean(obj.first_pass_peaks((end-3):end,:),1));
                                pks_ok = all(past_pk_errs<obj.pk_target_error);
                                % check widths are consistent
                                %width_errs = abs(obj.curr_widths-mean(obj.curr_widths));
                                %widths_ok =  all(width_errs<4);
                                widths_ok = true;
                                if pks_ok && widths_ok
                                    obj.found_candidate_peaks = true;
                                else
                                    obj.first_pass_peaks = [obj.first_pass_peaks ; obj.curr_locs];
                                    obj.windows = getWindows(obj.smp_freqs, obj.curr_locs, 10);
                                end
                            else
                                obj.first_pass_peaks = [obj.first_pass_peaks ; obj.curr_locs];
                                obj.windows = getWindows(obj.smp_freqs, obj.curr_locs, 10);
                                
                            end
                            
                        catch
                            obj.curr_locs = curr_guess(1:3:(end-1)); % if fit fails we still want peaks
                            obj.windows = getWindows(obj.smp_freqs, obj.curr_locs, obj.freq_window);
                            obj.found_candidate_peaks = false;
                        end
                    else
                        %obj.curr_locs = curr_guess(1:3:(end-1)); % can still use approximate peak locations to guess
                        %obj.windows = getWindows(obj.smp_freqs, obj.curr_locs, obj.freq_window);
                        obj.first_pass_peaks = [];
                        obj.found_candidate_peaks = false;
                    end
                else
                    obj.first_pass_peaks = [];
                    obj.found_candidate_peaks = false;
                end
                
            end
            
            if obj.found_candidate_peaks
                target_amp = 0.01;
                guessed_func = zeros(length(obj.smp_freqs),1);
                for i=1:length(obj.curr_locs)
                    guessed_func=guessed_func+target_amp*(0.5*pi*obj.curr_widths(i))*...
                        obj.basis_kernel(obj.smp_freqs,obj.curr_widths(i),...
                        obj.curr_locs(i));
                end
                [~, curr_num_pks, curr_guess] = getFitGuess(obj.smp_freqs,...
                    guessed_func, obj.contrast);
                if curr_num_pks == obj.target_num_pks
                    obj.num_second_pass = obj.num_second_pass+1;
                    try
                        [obj.fit, curr_params, ~, ~, curr_conf, obj.L_func] = ...
                            lorentzian_fit_mtx(obj.smp_freqs',...
                            (obj.projs-mean(obj.projs))', 2, 2, curr_num_pks, obj.smp_mtx,...
                            curr_guess);
                        [obj.curr_locs, obj.curr_locs_unc] = getFitVals(curr_params, ...
                            curr_conf, 'Peak');
                        [obj.curr_widths, obj.curr_widths_unc] = getFitVals(curr_params, ...
                            curr_conf, 'Width');
                        
                        % look at uncertainties and consistency in second pass peaks
                        obj.second_pass_peaks = [obj.second_pass_peaks ; obj.curr_locs];
                        if size(obj.second_pass_peaks,1)>=2
                            % compare second pass with first
                            past_pk_errs = abs(obj.second_pass_peaks(end,:)...
                                -mean(obj.second_pass_peaks((end-1):end,:),1));
                            pks_ok = all(past_pk_errs<obj.pk_target_error);
                            % check widths are consistent
                            width_errs = abs(obj.curr_widths-mean(obj.curr_widths));
                            widths_ok =  all(width_errs<5); %3);
                            %obj.num_second_pass = obj.num_second_pass+1;
                            if pks_ok && widths_ok
                                obj.fit_widths = true;
                                obj.fit_converged = true;
                                obj.fit = obj.curr_fit;
                                obj.pk_locs = obj.curr_locs;
                                obj.pk_unc = obj.curr_locs_unc;
                                obj.width = obj.curr_widths;
                                obj.width_unc = obj.curr_widths_unc;
                                
                            else
                                obj.second_pass_peaks = [obj.second_pass_peaks ; obj.curr_locs];
                                obj.windows = getWindows(obj.smp_freqs, obj.curr_locs, 10);
                            end
                        else
                            obj.second_pass_peaks = [obj.second_pass_peaks ; obj.curr_locs];
                            obj.windows = getWindows(obj.smp_freqs, obj.curr_locs, 10);
                            %obj.num_second_pass = obj.num_second_pass+1;
                        end
                    catch
                        obj.found_candidate_peaks = false;
                        obj.fit_converged = false;
                        obj.first_pass_peaks = [];
                        obj.opts.init = 1;
                    end
                else
                    obj.found_candidate_peaks = false;
                    obj.fit_converged = false;
                    obj.first_pass_peaks = [];
                    obj.opts.init = 1;
                    obj.fit_widths = false;
                end
            end
            
            
            %             if obj.fit_widths &&  (length(obj.curr_locs)==8)
            %                 %last pass
            %                 widths = zeros(obj.num_width_bins*obj.target_num_pks,1);
            %                 obj.basis_mtx = zeros(obj.num_pts, obj.num_width_bins*obj.target_num_pks);
            %                 for peak_idx=1:obj.target_num_pks
            %                     curr_indices = ((peak_idx-1)*obj.num_width_bins+1):(obj.num_width_bins*peak_idx);
            %                     widths(curr_indices)=linspace(obj.curr_widths(peak_idx)-0.2*obj.width_window, ...
            %                         obj.curr_widths(peak_idx)+obj.width_window,...
            %                         obj.num_width_bins);
            %                     obj.basis_mtx(:,curr_indices) = get_basis_mtx(obj.basis_kernel, obj.smp_freqs, obj.curr_locs(peak_idx), ...
            %                         widths(curr_indices));
            %                 end
            %                 obj.basis_widths = widths;
            %                 sensing_mtx = obj.smp_mtx*obj.basis_mtx;
            %
            %
            %
            %                 tau = 0.07*max(1e-4*max(abs(sensing_mtx'*(obj.projs-mean(obj.projs)))),obj.est_std*sqrt(log(length(obj.smp_freqs))));
            %                 delx_mode = 'mil';
            %                 %err_fun = @(z) (norm(x-z)/norm(x))^2;
            %                 %in = [];
            %                 obj.in.tau = tau;
            %                 obj.in.nonneg = 1; % Positivity constraint flag
            %                 obj.in.delx_mode = delx_mode;
            %                 obj.in.debias = 0;
            %                 obj.in.verbose = 0;
            %                 obj.in.plots = 0;
            %                 obj.in.record = 0;
            %                 %in.err_fun = err_fun;
            %                 x = l1homotopy(sensing_mtx,obj.projs-mean(obj.projs),obj.in);
            %                 obj.sparse_rep = x.x_out;
            %
            % %
            % %                 [obj.sparse_rep, ~] = TVAL3(sensing_mtx,obj.projs-mean(obj.projs),...
            % %                     size(obj.basis_mtx,2),1,obj.opts);
            %
            %
            %
            %
            %
            %
            % %                     test_prob_size = size(obj.basis_mtx,2);
            % %                     cvx_begin quiet
            % %                     cvx_precision best
            % %                     variable x(test_prob_size)
            % %                     dual variable y;
            % %                     dual variable z;
            % %                     minimize(norm(sensing_mtx * x-obj.projs+mean(obj.projs),1))
            % %                     %minimize(norm(x,1))
            % %                     subject to
            % %                     %y: norm((sensing_mtx*x-projections),2) <= (0.006)*sqrt(num_projs)*norm(projections,2);
            % %                     %y: sensing_mtx * x == obj.projs-mean(obj.projs);
            % %                     z: x >= 0;
            % %                     cvx_end
            % %                     obj.sparse_rep = x;
            %
            %
            %                 obj.reconstruction =  obj.basis_mtx*obj.sparse_rep;
            %
            %                 if obj.display
            %                     if sum(obj.fit)==0
            %                         plot(obj.ax, obj.smp_freqs, obj.reconstruction)
            %                     else
            %                         plot(obj.ax, obj.smp_freqs, obj.fit-mean(obj.fit))
            %                     end
            %                     title(num2str(obj.measurement_number))
            %                     drawnow
            %                 end
            %
            %                  [~, curr_num_pks, curr_guess] = getFitGuess(obj.smp_freqs,...
            %                     obj.reconstruction, obj.contrast);
            %                 if curr_num_pks == obj.target_num_pks
            %                     try
            %                         [obj.curr_fit, curr_params, ~, ~, curr_conf, obj.L_func] = ...
            %                             lorentzian_fit_mtx(obj.smp_freqs',...
            %                             obj.projs', 2, 2, curr_num_pks, obj.smp_mtx,...
            %                             curr_guess);
            %                         [obj.curr_locs, obj.curr_locs_unc] = getFitVals(curr_params, ...
            %                             curr_conf, 'Peak');
            %                         [obj.curr_widths, obj.curr_widths_unc] = getFitVals(curr_params, ...
            %                             curr_conf, 'Width');
            %                         obj.fit_converged = true;
            %
            %                         %obj.Continue = false;
            %                         obj.fit = obj.curr_fit;
            %                         obj.pk_locs = obj.curr_locs;
            %                         obj.pk_unc = obj.curr_locs_unc;
            %                         obj.width = obj.curr_widths;
            %                         obj.width_unc = obj.curr_widths_unc;
            %                         if obj.display
            %                             plot(obj.ax, obj.smp_freqs, obj.fit-mean(obj.fit))
            %                             title(num2str(obj.measurement_number))
            %                         end
            %
            %
            %
            %                     catch
            %                         obj.fit_converged = false;
            %                         obj.first_pass_peaks = [];
            %                         obj.opts.init = 1;
            %                         obj.fit_widths = false;
            %                     end
            %                 else
            %                     obj.fit_converged = false;
            %                     obj.first_pass_peaks = [];
            %                     obj.opts.init = 1;
            %                     obj.fit_widths = false;
            %                 end
            %             end
            
            
            if obj.get_movie && obj.flag_start_reconstructions
                obj.animation(obj.measurement_number-obj.num_initial_projs+1) = getframe(gcf);
                writeVideo(obj.video, obj.animation(obj.measurement_number-obj.num_initial_projs+1));
            end
            
            if obj.flag_start_reconstructions
                if sum(obj.fit)==0
                    obj.all_reconstructions = [obj.all_reconstructions; obj.reconstruction'];
                else
                    obj.all_reconstructions = [obj.all_reconstructions; obj.fit'];
                end
                obj.indice =  [obj.indice; obj.measurement_number];
                obj = obj.adaptive_method();
                if  obj.found_candidate_peaks && ~obj.save_statistics
                    obj.save_statistics = true;
                    obj.start_statistics = obj.measurement_number;
                end
            end
            
            if (obj.fit_converged) %&& obj.get_err
                obj.actual_err_locations = [obj.actual_err_locations; mean(abs(obj.pk_locs-obj.sim.peak_locs))];
                obj.actual_err = [obj.actual_err; sqrt( obj.actual_err_locations(end).^2 + mean(obj.curr_locs_unc).^2)];
                obj.meas_num = [obj.meas_num; obj.measurement_number];
                obj.all_err = [obj.all_err ; obj.actual_err(end)];

%                 obj.actual_err = [obj.actual_err; mean(abs(obj.pk_locs-obj.sim.peak_locs))];
%                 obj.meas_num = [obj.meas_num; obj.measurement_number];
            end
        end
        
        function obj = adaptive_method(obj)
%             num_peaks = length(obj.curr_locs);
%             prob_to_poll_peak = (num_peaks*obj.freq_window)/...
%                 (obj.smp_freqs(end)-obj.smp_freqs(1)); %one frequency
%             prob_to_poll_peak = 1-(1-prob_to_poll_peak)^obj.num_sample_freqs; %at least one frequency
%             
%             if ~obj.stop_polling_peaks && rand<=prob_to_poll_peak^(1/1.1)
%                 obj.poll_peaks_randomly = true;
%             else
%                 obj.poll_peaks_randomly = false;
%             end
%             
%             if ~obj.stop_polling_peaks && obj.found_candidate_peaks
%                 obj.poll_peaks_randomly = false;
%                 obj.stop_polling_peaks = true;
%             end
            
            obj.poll_peaks_randomly = false;
            if ~obj.poll_peaks_randomly
                prob_to_poll_window = (obj.window_end-obj.window_start)/...
                    (obj.smp_freqs(end)-obj.smp_freqs(1)); %one frequency
                prob_to_poll_window = 1-(1-prob_to_poll_window)^obj.num_sample_freqs;
                if rand <= (prob_to_poll_window)^(1/2)
                    obj.poll_window = true;
                else
                    obj.poll_window = false;
                end
            else
                obj.poll_window = false;
            end
            %obj.poll_window = false;
            
            
        end
    end
end


%% SUB FUNCTIONS

function y = lorentzian_kernel(varargin)
% Different choices of lorentzian kernel
% OUTPUT is an anonymous function whose form is a 1-D lorentzian with
% properties determined by varargin
% INPUT is a string:
% 'ZMN' or an empty string for Zero-Mean Normalized Lorentzian
% 'ZMU' for Zero-Mean Unormalized
% 'N' for Normalized
% 'U' for Unnormalized
switch length(varargin)
    case 0
        y = @(x,w,x0) (2/(pi*w)).*(1./(1+(2.*(x-x0)/w).^2))-...
            mean((2/(pi*w)).*(1./(1+(2.*(x-x0)/w).^2)));
        
    case 1
        if ~ischar(varargin{1})
            error('Choice for Lorentzian Kernel should be "ZMU", "ZMN", "U", or "N".');
        end
        switch varargin{1}
            case 'ZMU'
                y = @(x,w,x0) (1./(1+(2.*(x-x0)./w).^2))-...
                    mean((1./(1+(2.*(x-x0)./w).^2)));
            case 'ZMN'
                %                 y = @(x,w,x0) (2/(pi*w)).*(1./(1+(2.*(x-x0)./w).^2))-...
                %                     mean((2/(pi*w)).*(1./(1+(2.*(x-x0)./w).^2)));
                L_n = @(x,w,x0) (2/(pi*w)).*(1./(1+(2.*(x-x0)./w).^2));
                y = @(x,w,x0) L_n(x,w,x0)-(2./(pi*(max(x)-min(x)))).*atan((max(x)-min(x))./w);
                
            case 'U'
                y = @(x,w,x0) (1./(1+(2.*(x-x0)./w).^2));
            case 'N'
                y = @(x,w,x0) (2/(pi*w)).*(1./(1+(2.*(x-x0)./w).^2));
            otherwise
                error('Choice for Lorentzian Kernel should be "ZMU", "ZMN", "U", or "N".');
        end
    otherwise
        error('Choice for Lorentzian Kernel should be "ZMU", "ZMN", "U", or "N".');
end
end


function basis_mtx = get_basis_mtx(kernel, sample_pts, basis_pts, width_pts)
% Calculates the basis matrix based on 'kernel'
if length(width_pts)==1
    basis_mtx = zeros(length(sample_pts), length(basis_pts));
    for i=1:length(basis_pts)
        basis_mtx(:,i) = kernel(sample_pts, width_pts, basis_pts(i));
    end
else
    if length(basis_pts)==1
        basis_mtx = zeros(length(sample_pts), length(width_pts));
        for i=1:length(width_pts)
            basis_mtx(:,i) = kernel(sample_pts, width_pts(i), basis_pts);
        end
    else
        error('Either basis_pts or width_pts should be a scalar.')
    end
end
end

function basis_mtx = get_dual_basis_mtx(kernel, sample_pts, basis_pts, width_pts, mdpt)
% Calculates the basis matrix based on 'kernel'
if length(width_pts)==1
    basis_mtx = zeros(length(sample_pts), length(basis_pts));
    for i=1:length(basis_pts)
        basis_mtx(:,i) = kernel(sample_pts, width_pts, basis_pts(i))+...
            kernel(sample_pts, width_pts, 2*mdpt-(basis_pts(i)));
    end
else
    if length(basis_pts)==1
        basis_mtx = zeros(length(sample_pts), length(width_pts));
        for i=1:length(width_pts)
            basis_mtx(:,i) = kernel(sample_pts, width_pts(i), basis_pts);
        end
    else
        error('Either basis_pts or width_pts should be a scalar.')
    end
end
end







function [smooth_data, num_pks, initial_guess] = getFitGuess(t, data, ...
    contrast)
% Guess the parameters of the peaks in data based on variable 't'
smoothing = 0.01;
% num_pks = 0;
% while smoothing<0.08
%     smooth_data = smooth(data,smoothing,'loess'); %kills the hyperfine, we dont need it.
%     [~,locs,width,prominence] = findpeaks(smooth_data ,'MinPeakHeight', contrast);%in the last parameter insert the contrast of the smaiiest peak
%     num_pks=length(locs);
%     smoothing = smoothing+0.005;
%     if num_pks==8
%         break
%     end
% end
smooth_data = smooth(data,smoothing,'loess'); %kills the hyperfine, we dont need it.
[~,locs,width,prominence] = findpeaks(smooth_data ,'MinPeakHeight', contrast);%in the last parameter insert the contrast of the smaiiest peak
num_pks=length(locs);


if num_pks>8
    [prominence, idx] = sort(prominence, 'descend');
    prominence = prominence(1:8);
    width = width(idx(1:8));
    locs = locs(idx(1:8));
    num_pks=8;
end
initial_guess=zeros(1,num_pks+1);
for i=1:num_pks
    initial_guess(2*i-(2-i))=t(locs(i)); %peak location
    initial_guess(2*i-(2-i)+1)= width(i);     %width
    initial_guess(2*i-(2-i)+2)=-abs(smooth_data(locs(i))-min(smooth_data));%/abs(min(smooth_data));        %contrast
end
initial_guess(num_pks*3+1)=min(smooth_data);  %offset
end



function [vals, unc] = getFitVals(params, conf, parameter_name)
% After 'lorentzian_fit' get a value from the fit and its uncertainty.
switch parameter_name
    case 'Peak'
        idx=1;
    case 'Width'
        idx=2;
    case 'Amplitude'
        idx=3;
    case 'Offset'
        vals = params(end);
        unc = abs(conf(end,2)-conf(end,1));
        return
    otherwise
        error('Unknown Paramter of Lorentzian.')
end
if isempty(conf)
    unc=[];
else
    unc = abs(conf(idx:3:(end-1),2)-conf(idx:3:(end-1),1));
end

vals = params(idx:3:(end-1));
end

function [windows] = getWindows(freqs,peak_locs,pass_window)

windows = zeros(size(freqs));
for i=1:length(peak_locs)
    windows = windows+pulse(freqs,pass_window,peak_locs(i));
end
windows = double(windows>0);

end


function occurrences = getOccurs(mtx)
occurrences = sum(mtx~=0,1);
end
