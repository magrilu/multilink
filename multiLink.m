function [L,Dendro] = multiLink(X, P, modelStr, gricParam)
% This function provides a simple implementation of Multi-Link
% a multi-class multi-model fitting clustering algorithm based on
% preferences.
%
%
% INPUT:
%       * X:= input data points
%       * P:=  preference matrix
%       * modelStr:= string specifying the model(s) type used during
%       clustering
%       * gricParam:=  struct containing the GRIC parameters used for model selection
%         - gricParam.lambda1:= non negative wheight, the higher the more
%         model complexity is penalized
%         - gricParam.lambda2:= non negative wheight, the higher the more
%         model complexity is penalized
%         - gricParam.sigma:=
%
%
% OUTPUT:
%       * L:= Labeling vector representing the final clustering
%       * Dendro:= dendrogram matrix
%
%%
% If you use this code in your research please cite:
% Magri, Leveni, Boracchi, MultiLink: Multi-class Structure Recovery
% via Agglomerative Clustering and Model Selection. CVPR 2021
%
% Feel free to send suggestions and comments or to report bugs at luca.magri[at]polimi<dot>it
%%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
%IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
%DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
%TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
%OR THE USE OR OTHER DEALINGS IN THE SOFTWARE


%% for debugging/illustration purposes
do_debugPlot = 0; % display the evolution of the clustering
do_verbose   = 0;   % print information

%% debug options

if(do_debugPlot)
    cmap = brewermap(2*size(P,1),'Set2');
end
%% parsinfg
[isMergeableGric,cardmss] = parseMergeableGric(modelStr);
%  gric parameters
lambda1 = gricParam.lambda1;
lambda2 = gricParam.lambda2;
sigma   = gricParam.sigma;

%% build tanimoto distances between preferences
n  = size(P,1);
s0 = P*P'; % inner product
n0 = diag(s0); %norm
m0 = repmat(n0,1,n);
d0 = m0 + m0' - s0;
D  = 1 - s0./d0; % store the tanimoto distances between rows of P
%% single linkage based clustering
% checks on the distance matrix
if(size(D,1)~=size(D,2))
    % put the distance matrix in a squareform nxn
    D = squareform(D);
end
D = D + diag(inf.*ones(n,1)); % D distance between clusters

% preallocations
Dendro = zeros(n-1,3); % output dendrogram matrix.
L = 1:n;               % clusters labels
L = L(:);

%% main loop
s = 0;% iteration count
[dmin, i, j] = findMinDist(D,n);
while(dmin<Inf)
    
    cardi = sum(L==L(i));
    cardj = sum(L==L(j));
    % check if clusters of i and j can be merged
    if(cardi==1 && cardj==1)
        % if clusters are singleton, merge
        ok = true;
        isgric = false;
        isslnk = false;
    elseif(cardi< cardmss || cardj < cardmss)
        % if clusters aren't big enough to fit a model do Single Linkage
        % test
        ok = isMergeableSL(P, L, i, j);
        isgric = false;
        isslnk = true;
    else
        % perform gric test
        [ok, msScore, msOutput ] = isMergeableGric(X, L, i, j, lambda1 , lambda2, sigma);
        isgric = true;
        isslnk = false;
    end
    
    %% //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
    %                           display clustering evolution
    %  //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
    if(do_debugPlot && isgric)
        h = figure(99);
        clf;
        colorAccept ='g';
        hold all;
        
        if(isgric)
            if(ok)
                if(do_verbose)
                    disp(['merge by ', msScore.model])
                end
                title(['merge clusters by ', msScore.model],'fontsize',30,'fontname','Source Serif Pro')
                if(strcmp(msScore.model,'line'))
                    display_band( X,msOutput.mij(:) , sigma,colorAccept,0.7)
                elseif(strcmp(msScore.model,'circle'))
                    display_anulus( X, msOutput.mij(:), sigma,colorAccept,0.7)
                end
            else
                if(do_verbose)
                    disp(['merge rejected by ', msScore.model])
                end
                title(['\color{red}reject merge'],'fontsize',30,'fontname','Source Serif Pro')
            end         
        end
        bb = getBB(X,0.1);
        set(h,'color','white');
        displaySLMerge(X, L,i,j, n,cmap);
        xlim([bb.xmin,bb.xmax]);
        ylim([bb.ymin,bb.ymax]);
        drawnow;
        if(isslnk)
            if(ok)
                disp('SL merge accepted');
            else
                disp('SL merge rejected');
                pause(1)
            end
        end
    end
    % //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
    
    if(ok)
        if(do_verbose)
            fprintf('iteration: %i\n',s);
        end
        s = s+1;
        % update dendrogram
        Dendro(s,:) = [L(i) L(j) dmin];
        % update clustering
        L(L==L(i)) = n+s;
        L(L==L(j)) = n+s;
        
        % update cluster distance matrix
        D =  updateSlnkDist(D,i,j);
    else
        %fprintf('%i no merge\n',s);
        D(i,j) = Inf;
        D(j,i) = Inf;
    end
    % find minimum
    [dmin, i, j] = findMinDist(D,n);
    
end

Dendro(:,[1 2])=sort(Dendro(:,[1 2]),2);
L = grp2idx(L);
end




function [dmin, row, col] = findMinDist(D,n)
% find the minimum of a square distance matrix of size n
% the diagonal of the matrix should be set to Inf
% return the minimum value, its row and col subscript indices
[dmin,index] = min(D(:));
[row,col] = ind2sub([n,n],index);
end

function D =  updateSlnkDist(D,i,j)
% update the cluster distance matrix
% according to the single linkage rule
D(i,:) = min(D(i,:),D(j,:));
D(:,i) = D(i,:);
D(i,i) = Inf;
D(j,:) = Inf;
D(:,j) = Inf;

end