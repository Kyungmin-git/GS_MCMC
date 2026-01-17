function [sigma_k_log]= sigma_propose2(current_sigma_k_log)

global stepsize2
sigma_k_log=current_sigma_k_log+randn*stepsize2;

if sigma_k_log<-2
    sigma_k_log=-4-sigma_k_log;
elseif sigma_k_log>1
    sigma_k_log=2-sigma_k_log;
else
    
end

end
