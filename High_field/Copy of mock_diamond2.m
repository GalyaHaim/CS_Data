classdef mock_diamond2 < handle
    %MOCK_DIAMOND Creates Simulated NV Lineshape for Magnetic Sensing
    properties
        num_pts =450; % number of frequency points in simulation
        
       
        B_mag %= 50; % Gauss
        B_theta% =30; % Degrees
        B_phi %= 65; % Degrees
        
%         target= mock_diamond2;
 %   target.recalc_values
 %   target.sigma = 0.0004;  
% freqs = target.smp_freqs;
% y=target.getRaster();
% figure;plot(freqs,y)
        % Noise Simulation
        add_noise = true;%true;
        sigma = 0.00; %expressed as a fraction of the mean
        
        % Target Peak Amplitude - contrast
        peak_amp = 0.01%0.009;
        
        % Derived Quantities
        peak_locs = []; % will contain peak locations based on Magnetic Field
        smp_freqs = []; % where we measure in our simulation
        target = []; % lineshape we want to reconstruct
        peak_locs_from_fit = [];
        % model of what we actually measure
        sig = [];
        ref = [];
        
        % for combining 'sig' and 'ref' into a single array for compatibility
        % with adaptive_reconstruction class
        signal = zeros(1,1,2);
    end
    
    properties (Hidden)
        % Orientations of Diamond (in lab frame) assumed fixed
        diamond_unit_vecs = (1/sqrt(3))*[1 1 1;...
                                        1 -1 -1;...
                                        -1 1 -1;
                                        -1 -1 1];


%         XYZtoNV = 1/sqrt(3) * [1  1 -1 -1                               % transformation matrix from XYZ system to NV system
%                                1 -1  1 -1
%                                1 -1 -1  1];
% NVtoXYZ = 3/4 * XYZtoNV';    
        
        % Lorentzian Properties
        center_freq = 2877; % MHz
        width = 10; % MHz
        kernel = lorentzian_kernel('ZMN');
        move_center = false;
        center_offset = 0;
        shifts
        
        raster_flag = false;
        
        chk_noise = [];
    end
    
    methods
        function obj = mock_diamond2(varargin)
             obj.num_pts =650; % number of frequency points in simulation
        obj.B_mag = 110; % Gauss
        obj.B_theta =33; % Degrees
        obj.B_phi = 60; % Degrees
        
            
            if ~isempty(varargin)
                if length(varargin)==3
                    obj.raster_flag = true;
                    
                    obj.smp_freqs = varargin{1};
                    obj.num_pts = length(obj.smp_freqs);
                    
                    obj.sig = varargin{2};
                    obj.ref = varargin{3};
                    
                    obj.target = (obj.ref-obj.sig)/mean(obj.ref);                    
                    

                    
                    
                    [~, curr_num_pks, curr_guess] = getFitGuess(obj.smp_freqs,...
                        obj.target, 0.5*max(obj.target));
                    
                    if curr_num_pks ~=8
                       error('Too many peaks') 
                    end
                    [~, curr_params, ~, ~, curr_conf] = ...
                        lorentzian_fit(obj.smp_freqs',...
                        obj.target', 2, 2, curr_num_pks, ...
                        curr_guess);
                    
                    [obj.peak_locs, ~] = getFitVals(curr_params, ...
                        curr_conf, 'Peak');
                    
                    obj.peak_locs = obj.peak_locs';
                    
                    
                    
                    
                    
                elseif length(varargin)==4
                    obj.raster_flag = true;
                    obj.smp_freqs = varargin{1};
                    obj.sig = varargin{2};
                    obj.ref = varargin{3};
                    obj.sigma = varargin{4};
                    
                    obj.num_pts = length(obj.smp_freqs);
                    obj.target = (obj.ref-obj.sig);
                    
                    
                    obj.sig = varargin{2};
                    obj.ref = varargin{3};
                    
                    [~, curr_num_pks, curr_guess] = getFitGuess(obj.smp_freqs,...
                        obj.target, 0.5*max(obj.target));
                    
                    if curr_num_pks ~=8
                       error('Too many peaks') 
                    end
                    [~, curr_params, ~, ~, curr_conf] = ...
                        lorentzian_fit(obj.smp_freqs',...
                        obj.target', 2, 2, curr_num_pks, ...
                        curr_guess);
                    
                    [obj.peak_locs, ~] = getFitVals(curr_params, ...
                        curr_conf, 'Peak');
                    
                    obj.peak_locs = obj.peak_locs';

                    
                    
                    
                    
                else
                    error('Mock Diamond takes 0,3, or 4 args')
                    
                    
                end       
            else
                % mock_experiment with no arguments assumes default values
                % and creates a mock lineshape for running a Compressed
                % Sensing simulation
                
                % Calculate Peak Locations given the Magnetic Field
                B_vec = obj.pol2xyz();
                B_projs = sum(obj.diamond_unit_vecs.*B_vec,2);
                detunings = 2.8*B_projs; %distance from center in MHz, assume 2.8 MHz/G
                full_peak_window = max(obj.center_freq+detunings)-min(obj.center_freq-detunings);
                if obj.move_center
                    max_to_move = 20;
                    obj.center_offset = max_to_move*(2*rand-1);
                    %obj.center_offset =30;
                    obj.shifts = 0.5*(2*rand(8,1)-1);
                else
                    obj.shifts = 0;
                    obj.center_offset = 0;
                end
                
                obj.peak_locs = sort([obj.center_freq-detunings; ...
                    obj.center_freq+detunings])+obj.center_offset+...
                    obj.shifts;
                
                % Calculate Number of Points in Simulation
                % aim for a frequency spacing of ~2 MHz and that the peaks take
                % up about 60% of the full window
                fraction_peak_window = 0.6;
                df = 1; % MHz
                %full_peak_window = max(obj.peak_locs)-min(obj.peak_locs);
                if full_peak_window < 200
                    %failsafe if small magnetic field projections
                    full_peak_window = 200;
                end
                % amount to pad on either side
                pad_width = ((1-fraction_peak_window)/(2*fraction_peak_window))*full_peak_window;
                %obj.num_pts = ceil((full_peak_window+2*pad_width)/df);
                
                obj.smp_freqs = obj.center_freq+...
                    linspace(-(0.5*full_peak_window+pad_width), (0.5*full_peak_window+pad_width), obj.num_pts);
                % get actual lineshape and simulated signal and reference
                obj.target = obj.getLineShape(obj.smp_freqs);
                % get a mock raster scan
                %            [obj.sig, obj.ref] = obj.getRaster(obj.smp_freqs);
            end
            
            
            

        end
        
        function obj = recalc_values(obj)
            clear obj.peak_locs;
            B_vec = obj.pol2xyz();
                B_projs = sum(obj.diamond_unit_vecs.*B_vec,2);
                detunings = 2.8*B_projs; %distance from center in MHz, assume 2.8 MHz/G
                full_peak_window = max(obj.center_freq+detunings)-min(obj.center_freq-detunings);
                if obj.move_center
                    max_to_move = 20;
                    obj.center_offset = max_to_move*(2*rand-1);
                    %obj.center_offset =30;
                    obj.shifts = 0.5*(2*rand(8,1)-1);
                else
                    obj.shifts = 0;
                    obj.center_offset = 0;
                end
                
                obj.peak_locs = sort([obj.center_freq-detunings; ...
                    obj.center_freq+detunings])+obj.center_offset+...
                    obj.shifts;
                
                % Calculate Number of Points in Simulation
                % aim for a frequency spacing of ~2 MHz and that the peaks take
                % up about 60% of the full window
                fraction_peak_window = 0.6;
                df = 1; % MHz
                %full_peak_window = max(obj.peak_locs)-min(obj.peak_locs);
                if full_peak_window < 200
                    %failsafe if small magnetic field projections
                    full_peak_window = 200;
                end
                % amount to pad on either side
                pad_width = ((1-fraction_peak_window)/(2*fraction_peak_window))*full_peak_window;
                %obj.num_pts = ceil((full_peak_window+2*pad_width)/df);
                
                obj.smp_freqs = obj.center_freq+...
                    linspace(-(0.5*full_peak_window+pad_width), (0.5*full_peak_window+pad_width), obj.num_pts);
                % get actual lineshape and simulated signal and reference
                obj.target = obj.getLineShape(obj.smp_freqs);
                % get a mock raster scan
                %            [obj.sig, obj.ref] = obj.getRaster(obj.smp_freqs);
                
                
            
            y = obj.getRaster();
           [~, curr_num_pks, curr_guess] = getFitGuess(obj.smp_freqs,...
                        y, 0.001);
            %[yprime,params,resnorm,residual,conf]=lorentzian_fit(x,y,dhyp,Nhyp,Npeak,varargin)
if curr_num_pks ==8
     
    [fit, params,res,residual conf] = lorentzian_fit_lf(obj.smp_freqs,...  
                                y, 2, 2, 8, curr_guess);
           [curr_locs, curr_locs_unc] = getFitVals(params, conf, 'Peak');
            obj.peak_locs_from_fit = curr_locs;
else 
    disp('Problem with angles in mock_diamond')
end
          %[vals, unc] = getFitVals(params, conf, parameter_name)
        end
        
        function  vec = pol2xyz(obj, varargin)
            %pol2xyz Converts vector in spherical to rectangular
            %coordinates
            %   Assumes angles are given in degrees.
            if length(varargin)==2
                theta = varargin{1};
                phi = varargin{2};
                vec = [sind(theta)*cosd(phi) sind(theta)*sind(phi) cosd(theta)];
            elseif length(varargin)==3
                mag = varargin{1};
                theta = varargin{2};
                phi = varargin{3};
                vec = mag*[sind(theta)*cosd(phi) sind(theta)*sind(phi) cosd(theta)];
            elseif isempty(varargin)
                % defaults to magnetic field vector
                mag = obj.B_mag;
                theta = obj.B_theta;
                phi = obj.B_phi;
                vec = mag*[sind(theta)*cosd(phi) sind(theta)*sind(phi) cosd(theta)];
            else
                error('Vector should be given as [R,theta,phi] or [theta,phi].')
            end
        end
        
        function  lineshape = getLineShape(obj,freqs)
            %get the lineshape for the diamond at frequencies 'freqs'
            amp = obj.peak_amp/(2/(pi*obj.width));
            L = @(detuning) amp*obj.kernel(freqs,obj.width,detuning);
            lineshape = zeros(length(freqs),1);
            for i=1:length(obj.peak_locs)
                lineshape = lineshape + L(obj.peak_locs(i))';
            end
        end
        
        function  targ = getRaster(obj,varargin)
            %Get a Simulated Raster Scan with or without noise
            if isempty(varargin)
                freqs = obj.smp_freqs;
            else
                freqs = varargin{1};
            end
            
            
            if obj.raster_flag             
                targ = obj.target;
            else
                npts = length(freqs);
                raster_lineshape = obj.getLineShape(freqs);
                
                
                if obj.add_noise
                    sim_ref = ones(npts,1)+obj.sigma*randn(npts,1);
                    sim_sig = sim_ref-raster_lineshape+obj.sigma*randn(npts,1);
                else
                    sim_ref = ones(npts,1);
                    sim_sig = sim_ref-obj.target;
                end
                targ = (sim_ref-sim_sig)/mean(sim_ref);
            end
            
            
        end
        
        function  obj = getMeasurement(obj,sample_freqs)
            %Get a Simulated Raster Scan with or without noise
            num_samples = length(sample_freqs);
            indices = zeros(num_samples,1);
            for i=1:(num_samples)
                indices(i) = find(obj.smp_freqs==sample_freqs(i));
            end
            
            if obj.add_noise
                if obj.raster_flag
                    ref_noise = obj.sigma*randn*mean(obj.ref);
                    obj.signal(1,1,2) = mean(obj.ref(indices))+ref_noise;
                    obj.signal(1,1,1) = ref_noise+sum(obj.sig(indices))+obj.sigma*randn*mean(obj.ref);
                    obj.chk_noise = [obj.chk_noise; ref_noise];
                    
                else
                    obj.signal(1,1,2) = 1+obj.sigma*randn;
                    %obj.signal(1,1,1) = obj.signal(1,1,2)+obj.sigma*randn-mean(obj.target(indices));
                    obj.signal(1,1,1) = obj.signal(1,1,2)+obj.sigma*randn-sum(obj.target(indices));
                end
                
                
                
            else
                obj.signal(1,1,2) = 1;
                obj.signal(1,1,1) = 1-sum(obj.target(indices));
            end
            
        end
        
        function  [] = show_plots(obj)
            % plot
            figure
            hold on
            title('Orientation of Magnetic Field and Diamond')
            for i=1:length(obj.diamond_unit_vecs)
                x = [0 obj.diamond_unit_vecs(i,1)];
                y = [0 obj.diamond_unit_vecs(i,2)];
                z = [0 obj.diamond_unit_vecs(i,3)];
                line(x,y,z,'Color','Blue')
                scatter3(x(2),y(2),z(2),'b','filled')
            end
            B_vec = [sind(obj.B_theta)*cosd(obj.B_phi) sind(obj.B_theta)*sind(obj.B_phi) cosd(obj.B_theta)];
            quiver3(0,0,0,B_vec(1),B_vec(2),B_vec(3),'Color','Red')
            view(-28.2,10.2)
            
            figure
            plot(obj.smp_freqs,obj.target)
            grid on
            xlabel('Frequency (MHz)')
            ylabel('Zero-Mean Lineshape')
            title_str = ['Lineshape Model, ' num2str(obj.B_mag)...
                ' G, \theta=' num2str(obj.B_theta) char(176) ', \phi=' ...
                num2str(obj.B_phi) char(176)];
            title(title_str)
            xlim([obj.smp_freqs(1) obj.smp_freqs(end)])
            
        end
        
        
    end
end



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

