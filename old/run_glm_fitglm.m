%% GLM method using fitglm
%% Requires looping through each channel and itme point, very computation-heavy

% % Fit the GLM using IRLS
% tstart = tic;
% for iChan = 1%:nChan
%
%     disp('--------------------------------------------------------')
%     fprintf('          SUBJECT %g\n',iSub)
%     disp('--------------------------------------------------------')
%
%     parfor iFrame = 1:nFrames
%
%         % Extract the data for the current electrode and time point
%         Y = squeeze(data(iChan, iFrame, :));
%
%         h = lillietest(Y);
%         figure; histfit(Y);
%         if h, warning("Distribution is not normal. Use different distribution for GLM?"); end
%
%         % Debugging step: Check for NaN or Inf in the response variable
%         if any(isnan(Y)) || any(isinf(Y))
%             error('Response variable contains NaN or Inf at channel %d, time point %d', iChan, iFrame);
%         end
%
%         switch optimization
%
%             case 'IRLS'
%
%                 % Initialize IRLS w
%                 w = ones(nTrials, 1);
%
%                 for iter = 1:max_iter
%                     % Fit the GLM with current w
%                     mdl = fitglm(design_matrix_pca, Y, 'linear', 'Distribution', 'normal', 'w', w);
%
%                     % Compute residuals
%                     residuals_current = mdl.Residuals.Raw;
%
%                     % Update w (using Huber's weight function as an example)
%                     w_new = ones(size(residuals_current));
%                     abs_residuals = abs(residuals_current);
%                     w_new(abs_residuals > k) = k ./ abs_residuals(abs_residuals > k);
%
%                     % Check for convergence
%                     if max(abs(w_new - w)) < tol
%                         break;
%                     end
%
%                     % Update w for the next iteration
%                     w = w_new;
%                 end
%
%             case 'WLS'
%
%                 disp("Calculating weights using principal components projection...")
%                 [betas_limo, w, ~] = limo_WLS(design_matrix, Y);
%                 mdl = fitglm(design_matrix, Y, 'linear', 'Distribution', 'normal', 'w', w);
%         end
%
%         % Store the results
%         betas(iChan, iFrame, :) = mdl.Coefficients.Estimate;
%         residuals(iChan, iFrame, :) = mdl.Residuals.Raw;
%         R2(iChan, iFrame) = mdl.Rsquared.Ordinary;
%         % F(iChan, iFrame) = mdl.Coefficients.FStat(2:end); % Skip the intercept F-statistic
%         pvals(iChan, iFrame,:) = mdl.Coefficients.pValue(2:end); % Skip the intercept p-value
%     end
% end
% fprintf('Time to run GLM on one subject: %g (min)\n',round(toc(tstart)/60,1))
%
% % Save the results
% % save('GLM_IRLS_PCA_results.mat', 'betas', 'residuals', 'R2', 'F', 'pvals');

%% Same but with categorical variables instead of design matrix

% % Define the categorical design matrix
% condition_labels = cell(size(event_types)); % Create a cell array for categorical labels
% condition_labels(event_types == 2) = {'pleasant'};
% condition_labels(event_types == 4) = {'neutral'};
% condition_labels(event_types == 8) = {'unpleasant'};
% condition_labels = categorical(condition_labels);
%
% % Fit the GLM with current w
% tbl = table(condition_labels, Y, w, 'VariableNames', {'Condition', 'Response', 'w'});
% mdl = fitglm(tbl, 'Response ~ Condition', 'Distribution', 'normal', 'w', 'w');
%
% % Compute residuals
% residuals_current = mdl.Residuals.Raw;

