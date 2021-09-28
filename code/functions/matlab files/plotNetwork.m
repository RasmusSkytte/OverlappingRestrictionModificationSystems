function plotNetwork(ax, G, varargin)
    labels = [];
    if numel(varargin) > 0
        labels = varargin{1};
    end
    axes(ax);
    plot(G, 'NodeColor', 'r' ,'EdgeColor', 'k', 'Layout', 'Circle', 'LineWidth', sqrt(G.Edges.Weight ./ max(G.Edges.Weight)), 'MarkerSize', max(sqrt(G.degree), 1), 'NodeLabel', labels, 'NodeFontSize', 7)
end
