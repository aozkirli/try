function dist = likelihood(m,s,sigma)
%LIKELIHOOD Computes the likelihood distribution for stimulus 
%orientation s, given sensory noise sigma

% unnormalized likelihood
unnormalized_likelihood = @(x) exp( (-1/(2*sigma^2)) * (mod(s - x + 90,180)-90).^2);

% compute normalization constant by numerical integration
z = integral(unnormalized_likelihood,-90,90);

% normalize likelihood
dist = (1/z) * unnormalized_likelihood(m);

end