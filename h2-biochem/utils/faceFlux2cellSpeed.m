function speed = faceFlux2cellSpeed(G, faceFlux)
% Compute cell-centered speed from face fluxes without forming an AD matrix.

if ~isfield(G.cells, 'centroids')
    G = computeGeometry(G);
end

[cellNo, cellFaces] = getCellNoFaces(G);
neighbors = getNeighbourship(G, 'Topological', true);
signs = 2*(neighbors(cellFaces(:, 1), 1) == cellNo) - 1;
cellFlux = signs.*faceFlux(cellFaces);

offsets = G.faces.centroids(cellFaces, :) - G.cells.centroids(cellNo, :);
faceToCell = sparse(cellNo, 1:numel(cellNo), 1, ...
    G.cells.num, numel(cellNo));

speedSquared = 0;
for d = 1:size(offsets, 2)
    velocity = faceToCell*(cellFlux.*offsets(:, d));
    velocity = velocity./G.cells.volumes;
    speedSquared = speedSquared + velocity.^2;
end
speed = speedSquared.^0.5;
end
