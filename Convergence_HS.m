% TODO: studty the convergence of cl and cd by increasing nPanel. Just make a
% quick graph.

nPanelConvSequence = [4, 6, 8, 10:20:500];
cd_ConvSequence = zeros(size(nPanelConvSequence));
cl_ConvSequence = zeros(size(nPanelConvSequence));

for conv_idx = 1:length(nPanelConvSequence)
	nPanel = nPanelConvSequence(conv_idx);
	draft_Hess_Smith;
	cd_ConvSequence(conv_idx) = cd;
	cl_ConvSequence(conv_idx) = cl;
end