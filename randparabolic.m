function n = randparabolic(mean, envelope_width)
% Parabolic random number: returns a random number sampled from a parabolic
% probability distribution (concave down, with envelope_width defining the
% distance between the two points at which the parabola intersects the
% horizontal axis
if ~exist('mean','var') || isempty(mean)
    mean = 3;
end

if ~exist('envelope_width','var') || isempty(envelope_width)
    envelope_width = 4.23;
end

w = envelope_width;

% Probability Density Function
f = @(x) 3/(2*w)*(1 - 4/w^2*(x-mean).^2);
dx = w/(1e3); % step size
X = (mean-w/2):dx:(mean+w/2);

% Discretize f(x)
P = f(X);
P(X < mean - w/2) = 0;
P(X > mean + w/2) = 0;
P = 1/sum(P(:))*P; % normalize

% Cumulative sum
C = cumsum(P);
n = X(ceil(min(X) + sum(C(end)*rand > C)));
end