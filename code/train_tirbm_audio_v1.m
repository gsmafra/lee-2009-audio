function train_tirbm_audio_v1(ws, num_bases, spacing, pbias, pbias_lb, pbias_lambda, epsilon, l2reg, epsdecay, nPC, sigmaPC)

addpath ../voicebox

batch_size = 1;

% 1. Collect the spectrograms
tic
Pall = collect_speech_data2();
toc

% 2. Run PCA
[X startframe_list] = concatenate_speech_data(Pall, randsample(length(Pall), min(length(Pall), 1000)));
% X = X-repmat(mean(X,2), 1, size(X,2)); % don't subtract mean
Cov = X*X'/size(X,2);
[V D] = eig(Cov);

% 3. Try to reconstruct the original data with PCA components
numfeat = size(X,1);
E = V(:, numfeat:-1:numfeat-nPC+1);
S = diag(subvec(diag(D), numfeat:-1:numfeat-nPC+1));
Xpc = E'*X;
Xrec = E*Xpc;

Ewhiten = diag(sqrt(diag(S)+sigmaPC).^-1)*E';
Eunwhiten = E*diag(sqrt(diag(S)+sigmaPC));

Xrec2 = Eunwhiten*Ewhiten*X;

#figure, imagesc(X)
#figure, imagesc(Xrec2)
#break

% 4. Whitening (scaling)
Xw = Ewhiten*X;

% initialize the TIRBM parameters
sigma_start = 0.2;
sigma_stop = 0.2;

CD_mode = 'exp';
bias_mode = 'simple';

% Etc parameters
K_CD = 1;

% Initialization
W = [];
vbias_vec = [];
hbias_vec = [];
pars = [];

C_sigm = 1;

numchannels = size(Xw,1);

num_trials = 1;

% Initialize variables
if ~exist('pars', 'var') || isempty(pars)
    pars=[];
end

if ~isfield(pars, 'ws'), pars.ws = ws; end
if ~isfield(pars, 'num_bases'), pars.num_bases = num_bases; end
if ~isfield(pars, 'spacing'), pars.spacing = spacing; end

if ~isfield(pars, 'pbias'), pars.pbias = pbias; end
if ~isfield(pars, 'pbias_lb'), pars.pbias_lb = pbias_lb; end
if ~isfield(pars, 'pbias_lambda'), pars.pbias_lambda = pbias_lambda; end
if ~isfield(pars, 'bias_mode'), pars.bias_mode = bias_mode; end

if ~isfield(pars, 'std_gaussian'), pars.std_gaussian = sigma_start; end
if ~isfield(pars, 'sigma_start'), pars.sigma_start = sigma_start; end
if ~isfield(pars, 'sigma_stop'), pars.sigma_stop = sigma_stop; end

if ~isfield(pars, 'K_CD'), pars.K_CD = K_CD; end
if ~isfield(pars, 'CD_mode'), pars.CD_mode = CD_mode; end
if ~isfield(pars, 'C_sigm'), pars.C_sigm = C_sigm; end

if ~isfield(pars, 'num_trials'), pars.num_trials = num_trials; end
if ~isfield(pars, 'epsilon'), pars.epsilon = epsilon; end

if ~exist('W', 'var') || isempty(W)
    W = 0.01*randn(pars.ws, numchannels, pars.num_bases);
end

if ~exist('vbias_vec', 'var') || isempty(vbias_vec)
    vbias_vec = zeros(numchannels,1);
end

if ~exist('hbias_vec', 'var') || isempty(hbias_vec)
    hbias_vec = -0.1*ones(pars.num_bases,1);
end

fname_prefix = sprintf('../results/audio/tirbm_audio_LB_V1b_w%d_b%02d_pc%d_sigpc%g_p%g_pl%g_plambda%g_sp%d_%s_eps%g_epsdecay%g_l2reg%g_bs%02d_%s', ws, num_bases, nPC, sigmaPC, pbias, pbias_lb, pbias_lambda, spacing, CD_mode, epsilon, epsdecay, l2reg, batch_size, datestr(now, 30));

fname_save = sprintf('%s', fname_prefix);
fname_mat  = sprintf('%s.mat', fname_save);
mkdir(fileparts(fname_save));

initialmomentum  = 0.5;
finalmomentum    = 0.9;

error_history = [];
sparsity_history = [];

dWnorm_history_CD= [];
dWnorm_history_l2= [];

Winc=0;
vbiasinc=0;
hbiasinc=0;

for t=1:num_trials
	
    epsilon = pars.epsilon/(1+epsdecay*t);

    tic;
    ferr_current_iter = [];
    sparsity_curr_iter = [];

    for j=1:10
    
    	j
    	fflush(stdout);
    
        imdata = Ewhiten*Pall{ceil(rand()*length(Pall))};
        size(imdata)
       	size(W)
       	size(vbias_vec)
       	size(hbias_vec)
        return
        imdata = trim_audio_for_spacing_fixconv(imdata, ws, spacing);
        imdatatr = imdata';
        imdatatr = reshape(imdatatr, [size(imdatatr,1), 1, size(imdatatr,2)]);
    
        [ferr dW dh dv poshidprobs poshidstates negdata stat]= fobj_tirbm_CD_LB_sparse_audio(imdatatr, W, hbias_vec, vbias_vec, pars, CD_mode, bias_mode, spacing, l2reg);
        ferr_current_iter = [ferr_current_iter, ferr];
        sparsity_curr_iter = [sparsity_curr_iter, mean(poshidprobs(:))];

        if t<5,
            momentum = initialmomentum;
        else
            momentum = finalmomentum;
        end

        dWnorm_history_CD(j,t) = stat.dWnorm_CD;
        dWnorm_history_l2(j,t) = stat.dWnorm_l2;

        % update parameters
        Winc = momentum*Winc + epsilon*dW;
        W = W + Winc;

        vbiasinc = momentum*vbiasinc + epsilon*dv;
        vbias_vec = vbias_vec + vbiasinc;

        hbiasinc = momentum*hbiasinc + epsilon*dh;
        hbias_vec = hbias_vec + hbiasinc;
    end
    
    mean_err = mean(ferr_current_iter);
    mean_sparsity = mean(sparsity_curr_iter);

    if (pars.std_gaussian > pars.sigma_stop) % stop decaying after some point
        pars.std_gaussian = pars.std_gaussian*0.99;
    end
    error_history(t) = mean(ferr_current_iter);
    sparsity_history(t) = mean(sparsity_curr_iter);

    toc

    fprintf('epoch %d error = %g \tsparsity_hid = %g\n', t, mean(ferr_current_iter), mean(sparsity_curr_iter));
    save(fname_mat, 'W', 'pars', 't', 'vbias_vec', 'hbias_vec', 'error_history', 'sparsity_history', 'E', 'S', 'sigmaPC', 'Ewhiten', 'Eunwhiten', 'dWnorm_history_CD', 'dWnorm_history_l2');
    
    if 0 %mod(t, 10)==0
        fname_mat_timestamp  = sprintf('%s_%04dEPOCHS.mat', fname_save, t);
        save(fname_mat_timestamp, 'W', 'pars', 't', 'vbias_vec', 'hbias_vec', 'error_history', 'sparsity_history', 'E', 'S', 'sigmaPC', 'Ewhiten', 'Eunwhiten', 'dWnorm_history_CD', 'dWnorm_history_l2');
    end
    
    disp(sprintf('results saved as %s\n', fname_mat));

end
