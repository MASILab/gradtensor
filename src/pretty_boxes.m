function pretty_boxes()

	% Improve how boxplots look in MATLAB. When calling boxplot, use the
	% following Name-Value pairs for best results:
	% - 'Symbol', '.'
	% - 'OutlierSize', 12
	% Must be called directly after boxplot().
	%
	% Leon Cai
	% MASI Lab
	% March 28, 2021

	handles = get(get(gca, 'children'), 'children');   % Get the handles of all the objects
	tags = get(handles, 'tag');   % List the names of all the objects

	for i = 1:numel(handles) % Make lines bold
		    set(handles(i), 'LineWidth', 1);
	end

	boxes = 1:numel(tags); % Get indices of the box objects and draw patches
	boxes = boxes(strcmp(tags, 'Box'));
	for i = numel(boxes):-1:1 % iterate backwards through handles
		    patch(handles(boxes(i)).XData, handles(boxes(i)).YData, handles(boxes(i)).Color, 'FaceAlpha', 0.4, 'EdgeColor', 'none'); % make a patch
	end

end
