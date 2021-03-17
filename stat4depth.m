% statistics for difference  between ephys corrected depth vs. claimed
% depth during implanting
load('reviseLengthTip.mat');
load('insertLengthTip.mat');
revisedMask = logical(revisedLength(:,2:end));
diffs = revisedLength - insertLengthMat;
diffs = diffs(:,2:end);
revisedDiffs = diffs(revisedMask);
meanDiff = mean(revisedDiffs);
stdDiff = std(revisedDiffs);
% plot
histogram(revisedDiffs,50)
% make sure 95% confidence interval
x = norminv([0.05 0.95]);
xrange = meanDiff + stdDiff*x;
% [-440,510]
