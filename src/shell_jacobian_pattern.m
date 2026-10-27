function pattern = shell_jacobian_pattern(param)
%SHELL_JACOBIAN_PATTERN Conservative sparsity for all supported link models.
% Include the nonlocal rupture cascade and adhesion gate: every patch edge
% can affect damage throughout the patch and adhesion throughout the shell.
N = param.N;
Ne = param.Ne;
incidence = spones(param.edgeAccumI + param.edgeAccumJ);
nodeNeighbors = spones(incidence * incidence.' + speye(N));
nodeXYZ = kron(nodeNeighbors, ones(3));
nodeEdge = kron(incidence, ones(3,1));
nodeSelf = kron(speye(N), ones(3,1));
edgeSelf = speye(Ne);
patch = sparse(double(param.damageMask(:) > 0));
patchCoupling = patch * patch.';
damagePosition = spones(patchCoupling * nodeEdge.');
adhesionGate = sparse(ones(N,1)) * patch.';

Zpp = sparse(3*N,3*N);
Zpe = sparse(3*N,Ne);
Zpn = sparse(3*N,N);
Zep = sparse(Ne,3*N);
Zee = sparse(Ne,Ne);
Zen = sparse(Ne,N);
Zne = sparse(N,Ne);
pattern = spones([ ...
    Zpp, speye(3*N), Zpe, Zpe, Zpe, Zpn; ...
    nodeXYZ, nodeXYZ, nodeEdge, nodeEdge, nodeEdge, nodeSelf; ...
    nodeEdge.', nodeEdge.', edgeSelf, edgeSelf, edgeSelf, Zen; ...
    nodeEdge.', nodeEdge.', edgeSelf, edgeSelf, edgeSelf, Zen; ...
    damagePosition, Zep, Zee, Zee, patchCoupling + edgeSelf, Zen; ...
    nodeSelf.', nodeSelf.', Zne, Zne, adhesionGate, speye(N)]);
end
