function Y = computeResi( X, opts)
%SAMPLER function for sampling hypotheses for MCT - T-Linkage.
%
% opts.sampling indicates the policy of sampling:
%
% 'uniform' sampling, mss are randomly extracted.
%
% 'localized' sampling, euclidean distance is used to promote the
% extraction of closest mss. 'opts.quantile'is a parameter the controls
% probability of extracting close points. The smaller, the less the
% probability of extracting points more distant than quantile.
%
% 'nearest' as localized sampling promote the extraction of close mss.
% 'min_neigh' specifies the number of neighborood in which mss are
% extracted.
%
% opts.robust = 'x84' a robust refitting is performed once a new model
% hypotesis is instantiated. 'off' to turn off this feature.
%
% opts.geo = 1 if geometric distances are used to compute resiudals.
%
% opts.voting specifies the different voting functions used to express
% preferences of points.
%
%  opts.voting = 'binary' opts.epsi is used as cutoff to express binary
%  voting as in J-Linkage.
%
% opts.voting = 'exp' exponential voting
%
% opts.voting = 'gauss' gaussian voting, as exp but with 0 derivative in 0
% and a more gentle cutoff
%
% Please cite
% L. Magri, A. Fusiello; Fitting Multiple Heterogeneous Models by Multi-class Cascaded T-linkage CVPR, 2019.

if(~isfield(opts,'sampling'))
    opts.sampling = 'uniform';
    fprintf('\t uniform sampling is adopted\n');
end
if(~isfield(opts,'robust') )
    opts.robust = 'off';
    fprintf('\t robust refitting is turned off\n')
end
if(~isfield(opts,'geo'))
    opts.geo = 0;
end
if(~isfield(opts,'voting'))
    opts.voting = 'binary';
    fprintf('\t binary voting is adopted\n');
end

if(opts.geo==1)
    %fprintf('\t geometric distances are adopted\n');
end

model     = opts.model;
sampling  = opts.sampling;
m         = opts.m;

n = size(X,2);

switch sampling
    case 'uniform'
        sampling_id = 0;
    case 'localized'
        if( ~isfield(opts,'distance'))
            D = squareform(pdist(X','euclidean'));
        else
            D = opts.distance;
        end
        if( ~isfield(opts,'quantile'))
            qnt = 0.65;
        else
            qnt = opts.quantile;
        end
        sigma = quantile(D(:),qnt)/4;
        sampling_id = 1;
    case 'nearest'
        tree = KDTreeSearcher(X');
        
        num_neigh = min( opts.num_neigh, size(X,2));
        
        sampling_id = 2;
    otherwise
        sampling = 'uniform';
        sampling_id = 0;
        warning('unifrom sampling is adopted');
end




R = nan(n,m);


switch model
    
    
    case 'line'
        cardmss = 2;
        S = nan(cardmss,m);
        H = nan(3,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                % create model
                h1 = fit_line(mss);
                flg_degen = is_invalid_line(h1);
                cont_degen = cont_degen+1;
                
            end
            
            if(flg_degen == 1)
                inds = nan(1,cardmss);
                h1 = rand(3,1);
            end
            
            S(:,j) = inds;
            
            
            % compute residuals
            d1 = abs(h1(1).*X(1,:) + h1(2)*X(2,:) + h1(3));
            
            % store sampling result
            H(:,j) = h1;
            R(:,j) = d1;
        end
        
        
        
    case 'line_from_circle'
        center = opts.circle(1:2);
        % compute angular coefficient
        
        
        cardmss = 1;
        S = nan(cardmss,m);
        H = nan(3,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            %if(strcmp(sampling,'uniform'))
            if( sampling_id == 0)
                inds = randsample(n, cardmss,false);
                %elseif(strcmp(sampling,'nearest'))
            elseif(sampling_id == 2)
                seed = randsample(n, 1);
                bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                inds = datasample(bucket, cardmss,'Replace',false);
                %elseif(strcmp(sampling,'localized'))
            elseif(sampling_id == 1)
                inds = nan(1,cardmss);
                inds(1) = randsample(n, 1);
                w = exp(-(D(inds(1),:).^2)/sigma.^2);
                w(inds(1))=0;
                for i=2:cardmss
                    inds(i) = randsample(n,1,true,w);
                    w(inds(i)) = 0;
                end
            end
            S(:,j) = inds;
            mss = X(:,inds);
            % create model
            h1 = fit_line_from_circle( mss, center);
            
            % compute residuals
            %d1 = abs(h1(1).*X(1,:) + h1(2)*X(2,:) + h1(3));
            d1 = res_line(X,h1);
            H(:,j) = h1;
            R(:,j) = d1;
        end
        
        
    case 'circle'
        cardmss = 3;
        H = nan(3,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if(sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                flg_degen = iscolinear(mss(:,1),mss(:,2),mss(:,3));
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(3,1);
                h1 = rand(3,1);
            else
                h1 = fit_circle_taubin(mss);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = abs( sqrt(sum((X-repmat(h1(1:2),1,n)).^2,1))-h1(3));
            
            H(:,j) = h1;
            R(:,j) = d1;
            
        end
        
    case 'parabola'
        cardmss = 3;
        S = nan(cardmss,m);
        H = nan(3,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                % create model
                
                flg_degen = iscolinear(mss(:,1),mss(:,2),mss(:,3));
                cont_degen = cont_degen+1;
                %h1 = fit_parabola(mss);
                
                
            end
            
            if(flg_degen == 1)
                disp('nooo')
                inds = nan(3,1);
                h1 = rand(3,1);
            else
                h1 = fit_parabola(mss);
            end
            
            S(:,j) = inds;
            
            
            % compute residuals
            d1 = res_parabola(X,h1);
            
            
            
            % store sampling result
            H(:,j) = h1;
            R(:,j) = d1;
        end
        
        
    case 'fundamental'
        %disp('fundamental matrix fitting')
        
        %
        cardmss = 8;
        H = nan(9,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                f_tentative = fit_fm(mss);
                flg_degen = is_fundamental_degen( mss, f_tentative , cardmss );
                %                 if(flg_degen==1)
                %                     disp('degenerate')
                %                 end
                
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(cardmss,1);
                h1 = rand(9,1);
            else
                h1 = f_tentative;
            end
            S(:,j) = inds;
            
            % compute residuals
            if(opts.geo==1)
                d1 = res_fm_geo(X, h1);
            else
                d1 = res_fm(X, h1);
            end
            
            
            H(:,j) = h1;
            R(:,j) = d1;
            
            
        end
        
    case 'affine_fundamental'
        cardmss = 4;
        H = nan(9,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                f_tentative = fit_fm_affine(mss);
                flg_degen = 0; %is_homography_degen( mss, f_tentative , cardmss );
                %                 if(flg_degen==1)
                %                     disp('degenerate')
                %                 end
                
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(cardmss,1);
                h1 = rand(9,1);
            else
                h1 = f_tentative;
            end
            S(:,j) = inds;
            
            % compute residuals
            if(opts.geo==1)
                d1 = res_fm_geo(X, h1);
            else
                d1 = res_fm(X, h1);
            end
            % refit model
            
            
            H(:,j) = h1;
            R(:,j) = d1;
            
        end
        
        
        
        %% %%%%
    case 'homography'
        cardmss = 4;
        H = nan(9,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 50;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                if(~validateMSS_homography(X, inds))
                    flg_degen = 1 ;
                else
                    h_tentative = fit_homography(mss);
                    flg_degen = is_homography_degen(X, h_tentative, inds);
                end
                %if(flg_degen==1)
                %    disp('degenerate')
                %end
                
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                disp('rand')
                %inds = randsample(n, cardmss,false);
                % mss = X(:,inds);
                % h1 = fit_homography(mss);
                inds = nan(cardmss,1);
                h1 = rand(9,1);
            else
                h1 = h_tentative;
            end
            S(:,j) = inds;
            
            % compute residuals
            if(opts.geo == 1)
                d1 = res_homography_geo(X, h1);
            else
                d1 = res_homography(X,h1);
            end
            
            H(:,j) = h1;
            R(:,j) = d1;
            
        end
        
        
        
        %%
        
        
    case 'homography_from_fund'
        
        f = opts.fund;
        F = reshape(f,[3,3]);
        e2 = epipole(F');%null(F');
        A = star(e2)*F;
        
        cardmss = 3;
        H = nan(9,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 20;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                x1 = X(1:3,inds);
                x2 = X(4:6,inds);
                flg1 = iscolinear(x1(:,1),x1(:,2),x1(:,3),'h')|| iscolinear(x2(:,1),x2(:,2),x2(:,3),'h');
                if(~flg1)
                    
                    h_tentative = fit_homography_from_fund(mss,A,e2);
                    
                    flg_degen = is_homography_from_fund_degen(X, h_tentative, F, inds);
                end
                
                %if(flg_degen==1)
                %    disp('degenerate')
                %end
                
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(cardmss,1);
                h1 = rand(9,1);
            else
                
                h1 = h_tentative;
            end
            S(:,j) = inds;
            
            % compute residuals
            
            d1 = res_homography(X, h1);
            
            
            H(:,j) = h1;
            R(:,j) = d1;
            
        end
        
        
    case 'flat'
        
        delta = opts.dim_subspace;
        cardmss = delta+1;
        f =size(X,1);
        S = nan(cardmss,m);
        H = nan(f+f*f,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                elseif(sampling_id == 1)
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                % create model
                h1 = fit_flat(mss,delta);
                flg_degen = is_invalid_flat(mss,delta);
                cont_degen = cont_degen+1;
                
            end
            
            if(flg_degen == 1)
                disp('nooo')
                inds = nan(3,1);
                h1 = rand(3,1);
            end
            
            S(:,j) = inds;
            
            
            % compute residuals
            d1 = res_flat(X,h1);
            
            
            % store sampling result
            H(:,j) = flatStruct2Vect(h1);
            R(:,j) = d1;
            
        end
        
    case 'plane'
        H = nan(4,m);
        cardmss = 3;
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                if(strcmp(sampling,'uniform'))
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                elseif(strcmp(sampling,'nearest'))
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss);
                elseif(strcmp(sampling,'localized'))
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                flg_degen = iscolinear(mss(:,1),mss(:,2),mss(:,3));
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(3,1);
                h1 = rand(3,1);
            else
                h1 = fit_plane(mss);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = res_plane(X,h1);
            
            
            H(:,j) = h1;
            R(:,j) = d1;
            
            
        end
        
    case 'cylinder'
        disp('refit not supported')
        H = nan(7,m);
        cardmss = 2;
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                if(strcmp(sampling,'uniform'))
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                elseif(strcmp(sampling,'nearest'))
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss);
                elseif(strcmp(sampling,'localized'))
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                flg_degen = is_cylinder_degen(mss);
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(2,1);
                h1 = rand(7,1);
            else
                h1 = fit_pc_cylinder(mss);
                w0 = h1(4:6); % direction of the cylinder;
                h1 = convertToFiniteCylinder(h1,X);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = res_pc_cylinder(X,h1);
            
            
            H(:,j) = h1;
            R(:,j) = d1;
            
        end
    case 'cylinder_axis'
        disp('cylinder with fixed direction')
        w0 = opts.axis;
        H = nan(7,m);
        cardmss = 3;
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                if(strcmp(sampling,'uniform'))
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                elseif(strcmp(sampling,'nearest'))
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss);
                elseif(strcmp(sampling,'localized'))
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                flg_degen = is_cylinder_degen(mss);
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(2,1);
                h1 = rand(7,1);
            else
                h1 = fit_pc_cylinder_ls_axis(mss,w0);
                h1 = convertToFiniteCylinder(h1,X);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = res_pc_cylinder(X,h1);
            
            
            H(:,j) = h1;
            R(:,j) = d1;
            
        end
    case 'sphere'
        disp('refit not supported')
        H = nan(4,m);
        cardmss = 4;
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                if(strcmp(sampling,'uniform'))
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                elseif(strcmp(sampling,'nearest'))
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss);
                elseif(strcmp(sampling,'localized'))
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                flg_degen = abs(det([mss(1:3,:); ones(1,4)]))<1e-8;
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(3,1);
                h1 = rand(3,1);
            else
                h1 = fit_sphere(mss);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = res_sphere(X,h1);
            
            
            H(:,j) = h1;
            R(:,j) = d1;
            
            
        end
    case 'affinity'
        cardmss = 3;
        H = nan(9,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 50;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                if(~validateMSS_affine(X,inds))
                    flg_degen = 1 ;
                else
                    h_tentative = fit_affinity(mss);
                    flg_degen = validateTheta_homography(X, h_tentative, inds);
                end
                %if(flg_degen==1)
                %    disp('degenerate')
                %end
                
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                %disp('rand')
                %inds = randsample(n, cardmss,false);
                % mss = X(:,inds);
                % h1 = fit_homography(mss);
                inds = nan(cardmss,1);
                h1 = rand(9,1);
            else
                h1 = h_tentative;
            end
            S(:,j) = inds;
            
            % compute residuals
            if(opts.geo == 1)
                d1 = res_homography_geo(X, h1);
            else
                d1 = res_homography(X,h1);
            end
            
            H(:,j) = h1;
            R(:,j) = d1;
            
        end
    case 'lc'
        % sampling of lines
        optsLine.model = 'line';
        optsLine.sampling =  opts.sampling;
        optsLine.m = floor(m/2);
        optsLine.robust =  opts.robust;
        optsLine.voting = opts.voting';
        Yline = computeResi(X, optsLine);
        optsCircle.model = 'circle';
        optsCircle.sampling = opts.sampling;
        optsCircle.m = floor(m/2);
        optsCircle.robust = opts.robust;
        optsCircle.voting = opts.voting;
        Ycircle = computeResi(X, optsCircle);
        
        S = [Yline.S,Ycircle.S];
        R =  [Yline.R,Ycircle.R];
        H.line = Yline.H;
        H.circle = Ycircle.H;
    otherwise
        error('the model is not defined: sampler supports ''line'' and ''circle'', ''homography'' and ''fundamental''.')
end



Y.S = S;
Y.H = H;
Y.R = R;

end

