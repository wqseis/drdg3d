

% BC_FAULT = 100 for a single fault surface
% BC_FAULT = 100,101,102,... for multiple fault surfaces
BC_FAULT = 100;

% BC_FREE = 1 for a free surface
BC_FREE = 1;

FtoV = [1,3,2;1,2,4;2,3,4;1,4,3];% point outward

MAX_NUM_RECV = 100;