function Rot = rot_loc2glob(Rloc)
%ROT_LOC2GLOB Summary of this function goes here
%   Detailed explanation goes here

Rot = kron(speye(9),Rloc);
end

