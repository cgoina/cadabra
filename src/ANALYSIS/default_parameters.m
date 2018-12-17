function default = default_parameters

% MOST LIKELY FLY ORIENTATION AND PATH
default.correct_orient = 0; % correct fly orientation
default.correct_positions = 0; % find most likely fly path

% DEFAULTS for LUNGING
default.lunging.thresh = 0.7;
default.lunging.obj1.area_min = 0.1;
default.lunging.obj1.vel = [1 200];
default.lunging.obj2.vel = [0 45];
default.lunging.obj1.acc = [15 2000];
default.lunging.obj1.length = [0.7 5];
default.lunging.obj1.pos_change = [0.05 5];
default.lunging.der_distc = [-2.1 -0.15];
default.lunging.distc = [0.9 4];
% save('./defparams.par','default');
