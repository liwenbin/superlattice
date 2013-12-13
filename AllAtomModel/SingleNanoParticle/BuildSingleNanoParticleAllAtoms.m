%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUILDING THE LAMMPS DATA FILE OF GOLD NANOPARTICLE THIOL SYSTEM
%%% ALL ATOM MODEL, NAMELY EVERY ATOM OF ALKANETHIOL, INCLUDING
%%% CARBON, HYDROGEN AND SULFUR ARE ALL EXPLICITLY MODELLED
%%% WRITTEN BY WENBIN LI, MIT
%%% NOV. 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% The program output a LAMMPS data file 'data.np' for molecular
%%% dynamics simulation of alkanethiol self-assembly on individual
%%% gold nanocrystal surfaces.

%%% The program also output a CFG configuration files "gold_np_thiol.cfg,
%%% which can be visualized using Ju Li's Atomeye program:
%%% http://li.mit.edu/Archive/Graphics/A/

%%% Allinger et al's MM3 force field will be used to model alkanethiols
%%% Reference: Allinger, N. L., Yuh, Y. H., & Lii, J-H. (1989) 
%%% "Molecular Mechanics. The MM3 Force Field for Hydrocarbons. 1."
%%% J. Am. Chem. Soc. 111, 8551-8565. 

%%% lattice constant of gold
a=4.0782;

%%% supercell size
BoxSize = 100;

%%% supercell matrix
H = [BoxSize 0 0; 0 BoxSize 0; 0 0 BoxSize];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  number of carbon atoms per thiol molecules
%%%  here we work with octane-thiol (SC8H17)
%%%  hence natom_thiol_C = 8, but you can change it
%%%  for example if you work with dodecanethiol
%%%  then natom_thiol_C = 12
natom_thiol_C = 8;

natom_thiol_H = 2*natom_thiol_C+1;
natom_thiol_SC = natom_thiol_C + 1;
natom_thiol = natom_thiol_SC + natom_thiol_H;

%%%  put the thiol molecules initially on a sphere centered
%%%  on the gold nanocrystal. The sphere has radius R
R = input('radius of the sphere to put thiols (say = 30):');

%%%  the azimuthal angle and polar angle are divided evenly
%%%  forming a grid, the thiol molecules are then put on these grids
%%%  total number of thiol molecules will be the product
%%%  of nphi and ntheta
nphi   = 8;
ntheta = 17;

%%%% total number of thiols
n_thiol = nphi*ntheta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUILDING ICOSAHEDRAL GOLD NANOCRYSTAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  layer of gold atom in gold icosahedral nanocrystal
%%%  different number of gold layers will result in
%%%  different number of gold atoms in the nanocrystal
%%%  gold_layer = 6 results in a nanoparticle with 561 gold
%%%  atoms. The nanoparticle diameter will be around 3 nm

gold_layer = 6;
n_gold     = (10*(gold_layer^3) + 11*gold_layer)/3 -5*(gold_layer^2)-1;

%%%  matrix store gold atom coordinates
Pos_gold = zeros(n_gold, 3);

%%%  center of gold nanoparticle
center = [0.5 0.5 0.5];

%%%  parater used to build icosohedron
phi = (sqrt(5)+1)/2;

%%%  coordinates of the twelve vertices of icosahedron
V = zeros(12,3); 

V(1,:) = [1 0 phi];  V(2,:) = [-1 0 phi];  V(3,:) = [0 -phi 1];
V(4,:) = [phi -1 0]; V(5,:) = [phi 1 0];   V(6,:) = [0 phi 1];

for i = 7:12
  V(i,:) = -V(i-6,:);
end

%%% end points coordinates of the thirty eges, each edge has two end points
E = zeros(2, 3, 30);

E(1,:,1) = [1 0 phi];  E(2,:,1) = [-1 0 phi];
E(1,:,2) = [1 0 phi];  E(2,:,2) = [0 -phi 1];
E(1,:,3) = [1 0 phi];  E(2,:,3) = [phi -1 0];
E(1,:,4) = [1 0 phi];  E(2,:,4) = [phi 1 0];
E(1,:,5) = [1 0 phi];  E(2,:,5) = [0 phi 1];
E(1,:,6) = [0 phi 1];  E(2,:,6) = [phi 1 0];
E(1,:,7) = [0 phi 1];  E(2,:,7) = [0 phi -1];
E(1,:,8) = [0 phi 1];  E(2,:,8) = [-phi 1 0];
E(1,:,9) = [0 phi 1];  E(2,:,9) = [-1 0 phi];
E(1,:,10)= [phi 1 0];  E(2,:,10)= [phi -1 0];
E(1,:,11)= [phi 1 0];  E(2,:,11)= [1 0 -phi];
E(1,:,12)= [phi 1 0];  E(2,:,12)= [0 phi -1];
E(1,:,13)= [phi -1 0]; E(2,:,13)= [1 0 -phi];
E(1,:,14)= [1 0 -phi]; E(2,:,14)= [0 phi -1];
E(1,:,15)= [0 phi -1]; E(2,:,15)= [-phi 1 0];

for i = 16:30
  E(:,:,i) = -E(:,:,i-15);
end

%%% coordinates of the vertices of twenty faces

F = zeros(3, 3, 20);

F(1,:,1)  = [1 0 phi];  F(2,:,1)  = [-1 0 phi];  F(3,:,1)  = [0 -phi 1];
F(1,:,2)  = [1 0 phi];  F(2,:,2)  = [0 -phi 1];  F(3,:,2)  = [phi -1 0];
F(1,:,3)  = [1 0 phi];  F(2,:,3)  = [phi -1 0];  F(3,:,3)  = [phi 1 0];
F(1,:,4)  = [1 0 phi];  F(2,:,4)  = [phi 1 0];   F(3,:,4)  = [0 phi 1];
F(1,:,5)  = [1 0 phi];  F(2,:,5)  = [0 phi 1];   F(3,:,5)  = [-1 0 phi];
F(1,:,6)  = [0 phi 1];  F(2,:,6)  = [phi 1 0];   F(3,:,6)  = [0 phi -1];
F(1,:,7)  = [0 phi 1];  F(2,:,7)  = [0 phi -1];  F(3,:,7)  = [-phi 1 0];
F(1,:,8)  = [0 phi 1];  F(2,:,8)  = [-phi 1 0];  F(3,:,8)  = [-1 0 phi];
F(1,:,9)  = [phi 1 0];  F(2,:,9)  = [phi -1 0];  F(3,:,9)  = [1 0 -phi];
F(1,:,10) = [phi 1 0];  F(2,:,10) = [1 0 -phi];  F(3,:,10) = [0 phi -1];

for i = 11:20
  F(:,:,i) = -F(:,:,i-10);
end

%%%  absolute coordinates of the inner-most gold layer
scaling = a/sqrt(2*(1+phi*phi));

V = V*scaling;
E = E*scaling;
F = F*scaling;

%%%  coordinates of the centeral atom
Pos_gold(1,:) = center;

%%%  coordinates of the second innermost layer atoms

for i=1:12
   Pos_gold(i+1,1) = V(i,1)/BoxSize + center(1);
   Pos_gold(i+1,2) = V(i,2)/BoxSize + center(2);
   Pos_gold(i+1,3) = V(i,3)/BoxSize + center(3);
end

%%%  coordinates of the third innermost layer atoms
for i= 1:12
  Pos_gold(i+13,1) = V(i,1)*2/BoxSize + center(1);
  Pos_gold(i+13,2) = V(i,2)*2/BoxSize + center(2);
  Pos_gold(i+13,3) = V(i,3)*2/BoxSize + center(3);
end

%%%  coordinates of the fourth innermost layer atoms
for i= 1:30
  Pos_gold(i+25,1) = (E(1,1,i)+E(2,1,i))/BoxSize+center(1);
  Pos_gold(i+25,2) = (E(1,2,i)+E(2,2,i))/BoxSize+center(2);
  Pos_gold(i+25,3) = (E(1,3,i)+E(2,3,i))/BoxSize+center(3);
end

%%%  the other layers of atoms, which have vertice atoms, edges atoms and face atoms

for N = 4:gold_layer
  %%%  total atom up to N-th layer
  NatomToNow = (10*((N-1)^3) + 11*(N-1))/3 -5*((N-1)^2)-1;
  
  %%%  coordinates of the twelve vertices of the N-th layer 
  for i=1:12
    Pos_gold(NatomToNow+i,1) = V(i,1)*(N-1)/BoxSize + center(1);
    Pos_gold(NatomToNow+i,2) = V(i,2)*(N-1)/BoxSize + center(2);
    Pos_gold(NatomToNow+i,3) = V(i,3)*(N-1)/BoxSize + center(3);
  end
  
  NatomToNow = NatomToNow+12;
  
  %%%  coordinates of the edge atoms
  for i=1:30
    for j = 1: (N-2)
      Pos_gold(NatomToNow+(i-1)*(N-2)+j,1) = (E(1,1,i)*(N-1)+ (E(2,1,i)-E(1,1,i))*j)/BoxSize + center(1);
      Pos_gold(NatomToNow+(i-1)*(N-2)+j,2) = (E(1,2,i)*(N-1)+ (E(2,2,i)-E(1,2,i))*j)/BoxSize + center(2);
      Pos_gold(NatomToNow+(i-1)*(N-2)+j,3) = (E(1,3,i)*(N-1)+ (E(2,3,i)-E(1,3,i))*j)/BoxSize + center(3);
    end
  end
  
  NatomToNow = NatomToNow + 30*(N-2);
  
  %%%  coordinates of the face atoms
  for i = 1:20
    for m = 1:(N-3)
      for n = 1:(N-m-2)
        Pos_gold(NatomToNow+ (2*N-m-4)*(m-1)/2+ n, 1) = (F(1,1,i)*(N-1) + (F(2,1,i)*(N-1)-F(1,1,i)*(N-1))*m/(N-1) + ...
                                                         (F(3,1,i)*(N-1)-F(1,1,i)*(N-1))*n/(N-1))/BoxSize+center(1);
        Pos_gold(NatomToNow+ (2*N-m-4)*(m-1)/2+ n, 2) = (F(1,2,i)*(N-1) + (F(2,2,i)*(N-1)-F(1,2,i)*(N-1))*m/(N-1) + ...
                                                         (F(3,2,i)*(N-1)-F(1,2,i)*(N-1))*n/(N-1))/BoxSize+center(2);
        Pos_gold(NatomToNow+ (2*N-m-4)*(m-1)/2+ n, 3) = (F(1,3,i)*(N-1) + (F(2,3,i)*(N-1)-F(1,3,i)*(N-1))*m/(N-1) + ...
                                                         (F(3,3,i)*(N-1)-F(1,3,i)*(N-1))*n/(N-1))/BoxSize+center(3);
      end
    end
    NatomToNow = NatomToNow + (N-2)*(N-3)/2;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUIDING ICOSAHEDRON GOLD NANOPARTICLES FINISHED
%%% NEXT BUILD THIOL MOLECULES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  coordinates of the thiol molecules
Pos_chain = zeros(natom_thiol,3,n_thiol);

%%%  matrix for conversion from absolute to reduced coordinates
H1 = [1/BoxSize 0 0; 0 1/BoxSize 0; 0 0 1/BoxSize];

%%%  the following vectors are used for determing the coordinates of thiol molecules
ShiftVector1    = [1.005 0 1.51];
ShiftVector2    = [-0.919 0 1.224];
TranslateVector = [-0.071 0 2.49];

ShiftVector_H1  = [1.634  -0.89 1.528];
ShiftVector_H2  = [1.634   0.89 1.518];
ShiftVector_H3  = [-0.543 -0.89 2.716];
ShiftVector_H4  = [-0.543  0.89 2.716];
ShiftVector_CH3 = [0.604    0   0.908];

n_unit          = natom_thiol_C;

k = 1;

%%%  increment of azimuthal and polar angles
dphi = 2*pi/nphi; dtheta = pi/(ntheta+3);

for phi = dphi:dphi:2*pi
  for theta = 2*dtheta:dtheta:dtheta*(ntheta+1)
    Pos_chain(1,1,k) = BoxSize/2 + R*sin(theta)*cos(phi);
    Pos_chain(1,2,k) = BoxSize/2 + R*sin(theta)*sin(phi);
    Pos_chain(1,3,k) = BoxSize/2 + R*cos(theta);
    
    %%% rotation matrix    
    M = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)
         -sin(phi)         cos(phi)            0
         sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];

    %%% rotate the thiols
    Pos_chain(2,:,k) = Pos_chain(1,:,k) + ShiftVector1*M;
    Pos_chain(3,:,k) = Pos_chain(1,:,k) + (ShiftVector1+ShiftVector2)*M;
    Pos_chain(natom_thiol_SC+1,:,k) = Pos_chain(1,:,k) + ShiftVector_H1*M;
    Pos_chain(natom_thiol_SC+2,:,k) = Pos_chain(1,:,k) + ShiftVector_H2*M;
    Pos_chain(natom_thiol_SC+3,:,k) = Pos_chain(1,:,k) + ShiftVector_H3*M;
    Pos_chain(natom_thiol_SC+4,:,k) = Pos_chain(1,:,k) + ShiftVector_H4*M;
    
    for n = 1:(n_unit-2)/2
      Pos_chain(3+2*n-1,:,k) = Pos_chain(1,:,k) + (ShiftVector1+TranslateVector*n)*M;
      Pos_chain(3+2*n,:,k)   = Pos_chain(1,:,k) + (ShiftVector1+ShiftVector2+TranslateVector*n)*M;
      
      Pos_chain(natom_thiol_SC+ 1 + 4*n,:,k) = Pos_chain(natom_thiol_SC+1,:,k) + (TranslateVector*n)*M;
      Pos_chain(natom_thiol_SC+ 2 + 4*n,:,k) = Pos_chain(natom_thiol_SC+2,:,k) + (TranslateVector*n)*M;
      Pos_chain(natom_thiol_SC+ 3 + 4*n,:,k) = Pos_chain(natom_thiol_SC+3,:,k) + (TranslateVector*n)*M;
      Pos_chain(natom_thiol_SC+ 4 + 4*n,:,k) = Pos_chain(natom_thiol_SC+4,:,k) + (TranslateVector*n)*M;
    end
    
    Pos_chain(natom_thiol,:,k) = Pos_chain(natom_thiol_SC,:,k) + ShiftVector_CH3*M;
    
    %%% Conversion from absolute coordinates to reduced coordinates
    Pos_chain(:,:,k) = Pos_chain(:,:,k)*H1;
    
    k = k+1;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% COORDINATES DETERMINATION ENDS HERE, NEXT OUTPUT CFG COORDINATION FILE
%%% THE CFG COORDINATION FILE CAN BE VISUALIZED BY JU LI'S ATOMEYE PROGRAM
%%% http://li.mit.edu/Archive/Graphics/A/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atom_name_Au = 'Au';
atom_mass_Au = 196.967;

atom_name_S = 'S';
atom_mass_S = 31.972;

atom_name_C = 'C';
atom_mass_C = 12.000;

atom_name_H = 'H';
atom_mass_H = 1.008;

natom_total = n_gold + natom_thiol*n_thiol;

cfg_name = 'gold_np_thiol.cfg';
cfg = fopen(cfg_name, 'w');
 
fprintf(cfg, 'Number of particles = %d\n', natom_total);
fprintf(cfg, 'H0(1,1)= %.5g A\n', H(1,1));
fprintf(cfg, 'H0(1,2)= %.5g A\n', H(1,2));
fprintf(cfg, 'H0(1,3)= %.5g A\n', H(1,3));
fprintf(cfg, 'H0(2,1)= %.5g A\n', H(2,1));
fprintf(cfg, 'H0(2,2)= %.5g A\n', H(2,2));
fprintf(cfg, 'H0(2,3)= %.5g A\n', H(2,3));
fprintf(cfg, 'H0(3,1)= %.5g A\n', H(3,1));
fprintf(cfg, 'H0(3,2)= %.5g A\n', H(3,2));
fprintf(cfg, 'H0(3,3)= %.5g A\n', H(3,3));
fprintf(cfg, '.NO_VELOCITY.\n');
fprintf(cfg, 'entry_count = 3\n');

fprintf(cfg, '%g\n', atom_mass_Au);
fprintf(cfg, '%2s\n', atom_name_Au);

for i=1: n_gold
  fprintf(cfg, '%.5g %.5g %.5g\n', Pos_gold(i,1), Pos_gold(i,2), Pos_gold(i,3));
end

fprintf(cfg, '%g\n', atom_mass_S);
fprintf(cfg, '%2s\n', atom_name_S);

for k = 1:n_thiol
  fprintf(cfg, '%.5g %.5g %.5g\n', Pos_chain(1,1,k), Pos_chain(1,2,k), Pos_chain(1,3,k));
end
  
fprintf(cfg, '%g\n', atom_mass_C);
fprintf(cfg, '%2s\n', atom_name_C);

for k = 1:n_thiol
  for i = 2:natom_thiol_SC
    fprintf(cfg, '%.5g %.5g %.5g\n', Pos_chain(i,1,k), Pos_chain(i,2,k), Pos_chain(i,3,k));
  end
end

fprintf(cfg, '%g\n', atom_mass_H);
fprintf(cfg, '%2s\n', atom_name_H);

for k = 1:n_thiol
  for i = (natom_thiol_SC+1):natom_thiol
    fprintf(cfg, '%.5g %.5g %.5g\n', Pos_chain(i,1,k), Pos_chain(i,2,k), Pos_chain(i,3,k));
  end
end
  
fclose(cfg);

%convert to absolute coordinates
Pos_gold = Pos_gold*H;

for k = 1:n_thiol
  Pos_chain(:,:,k)=Pos_chain(:,:,k)*H;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  NEXT OUTPUT LAMMPS DATA FILE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  total number of atoms
total_atoms = n_gold + natom_thiol*n_thiol;

%%%  number of bonds per thiol
n_bonds_per_thiol = natom_thiol_C + natom_thiol_C*2 + 1;

%%%  total number of bonds
total_bonds       = n_bonds_per_thiol*n_thiol;

%%%  total number of angles
n_angles_per_thiol = natom_thiol_C*6;
total_angles       = n_angles_per_thiol*n_thiol;

%%%  total number of dihedrals
n_dihedrals_per_thiol = (natom_thiol_C-1)*9;
total_dihedrals       = n_dihedrals_per_thiol*n_thiol;

%%%  atom types:
%%%  Au-1
%%%  S-2
%%%  C-3
%%%  H-4
atom_types = 4;

%%%  bond types:
%%%  1  S-C
%%%  2  C-C
%%%  3  C-H
bond_types = 3;

%%%  angle types:
%%%  1  C-C-C
%%%  2  C-C-H  the middle C connects with two hydrogens
%%%  3  C-C-H  the middle C connects with three hydrogens
%%%  4  H-C-H  the middle C connects with two hydrogens
%%%  5  H-C-H  the middle C connects with three hydrogens
%%%  6  S-C-C
%%%  7  S-C-H
angle_types = 7;

%%%  dihedral types
%%%  1  C-C-C-C
%%%  2  CH2-CH2-CH2-H
%%%  3  H-CH2-CH2-H
%%%  4  S-C-C-C
%%%  5  S-C-C-H
%%%  6  CH2-CH2-CH3-H
%%%  7  H-CH2-CH3-H
dihedral_types = 7;

atom_mass_Au = 196.967;
atom_mass_S =  31.972;
atom_mass_C =  12.000;
atom_mass_H =  1.008;

%%%  file name
file_name = 'data.np';

fp = fopen(file_name, 'w');

fprintf(fp, 'Gold Nanoparticle - Thiol System\n');
fprintf(fp, '\n');
fprintf(fp, '        %d   atoms\n', total_atoms);
fprintf(fp, '        %d   bonds\n', total_bonds);
fprintf(fp, '        %d   angles\n', total_angles);
fprintf(fp, '        %d   dihedrals\n', total_dihedrals);
fprintf(fp, '\n');
fprintf(fp, '        %d   atom types\n', atom_types);
fprintf(fp, '        %d   bond types\n', bond_types);
fprintf(fp, '        %d   angle types\n', angle_types);
fprintf(fp, '        %d   dihedral types\n', dihedral_types);

fprintf(fp, '\n');
fprintf(fp, '  0.00  %.5g        xlo xhi\n', H(1,1));
fprintf(fp, '  0.00  %.5g        ylo yhi\n', H(2,2));
fprintf(fp, '  0.00  %.5g        zlo zhi\n', H(3,3));

fprintf(fp, '\n');
fprintf(fp, 'Masses\n');
fprintf(fp, '\n');

fprintf(fp, ' 1 %g\n', atom_mass_Au);
fprintf(fp, ' 2 %g\n', atom_mass_S);
fprintf(fp, ' 3 %g\n', atom_mass_C);
fprintf(fp, ' 4 %g\n', atom_mass_H);

%%%%%%%%%%%%%%%%%%%%%
%%%  PRINT ATOMS  %%%
%%%%%%%%%%%%%%%%%%%%%
 
fprintf(fp, '\n');
fprintf(fp, 'Atoms\n');
fprintf(fp, '\n');

%%% PRINT FORMAT: ATOM_ID MOLECULE_ID ATOM_TYPE POSITION_X POSITON_Y POSITION_Z

%%% print gold atoms
for i = 1: n_gold
  fprintf(fp, '   %d     1   1   %.5g  %.5g  %.5g\n', i, Pos_gold(i,1), Pos_gold(i,2), Pos_gold(i,3));
end

for k = 1:n_thiol
  %%% print sulfur atoms
  fprintf(fp, '   %d     %d   2  %.5g  %.5g  %.5g\n', n_gold+(k-1)*natom_thiol+1, 1+k, Pos_chain(1,1,k), Pos_chain(1,2,k), Pos_chain(1,3,k));
  %%% print carbon atoms
  for i=2:natom_thiol_SC
    fprintf(fp, '   %d     %d   3  %.5g  %.5g  %.5g\n', n_gold+(k-1)*natom_thiol+i, 1+k, Pos_chain(i,1,k), Pos_chain(i,2,k), Pos_chain(i,3,k));
  end
  %%% print hydrogen atoms
  for i=(natom_thiol_SC+1):natom_thiol
    fprintf(fp, '   %d     %d   4  %.5g  %.5g  %.5g\n', n_gold+ (k-1)*natom_thiol+i, 1+k, Pos_chain(i,1,k), Pos_chain(i,2,k), Pos_chain(i,3,k));
  end
end

%%%%%%%%%%%%%%%%%%%%%
%%%  PRINT BONDS  %%%
%%%%%%%%%%%%%%%%%%%%%

fprintf(fp, '\n');
fprintf(fp, 'Bonds\n');
fprintf(fp, '\n');

%%% PRINT FORMAT: BOND_NUMBER BOND_TYPE ATOM_ID1 ATOM_ID2

for k = 1: n_thiol
  
  %%% print S-C bonds
  fprintf(fp, '     %d   %d     %d     %d\n',  (k-1)*n_bonds_per_thiol+1, 1, ...
            n_gold+(k-1)*natom_thiol+1, n_gold+(k-1)*natom_thiol+2);
  
  %%% print C-C bonds
  for i = 2:natom_thiol_C
    fprintf(fp, '     %d   %d     %d     %d\n',  (k-1)*n_bonds_per_thiol+i, 2, ...
            n_gold+(k-1)*natom_thiol+i, n_gold+(k-1)*natom_thiol+i+1);
  end
  
  %%% print C-H bonds
  for i = 1:natom_thiol_C
    fprintf(fp, '     %d   %d     %d     %d\n',  (k-1)*n_bonds_per_thiol + natom_thiol_C + 2*i-1, 3, ...
            n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i-1);
    fprintf(fp, '     %d   %d     %d     %d\n',  (k-1)*n_bonds_per_thiol + natom_thiol_C + 2*i,   3,...
            n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i);
  end
  
  %%% print the last C-H bond in each thiol molecule
  fprintf(fp,   '     %d   %d     %d     %d\n',  (k-1)*n_bonds_per_thiol+ 3*natom_thiol_C + 1, 3, ...
          n_gold+(k-1)*natom_thiol + natom_thiol_C + 1,  n_gold + k*natom_thiol);
end


%%%%%%%%%%%%%%%%%%%%%%
%%%  Print Angles  %%%
%%%%%%%%%%%%%%%%%%%%%%

fprintf(fp, '\n');
fprintf(fp, 'Angles\n');
fprintf(fp, '\n');

%%% PRINT FORMAT: ANGLE_NUMBER ANGLE_TYPE ATOM_ID1 ATOM_ID2 ATOM_ID3

for k = 1: n_thiol
  
    % S-C-C
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 1, 6, ...
            n_gold+(k-1)*natom_thiol+1, n_gold+(k-1)*natom_thiol+2, n_gold+(k-1)*natom_thiol+3);
    % S-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 2, 7, ...
            n_gold+(k-1)*natom_thiol+1, n_gold+(k-1)*natom_thiol+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 1);
    % S-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 3, 7, ...
            n_gold+(k-1)*natom_thiol+1, n_gold+(k-1)*natom_thiol+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2);
    % C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 4, 2, ...
            n_gold+(k-1)*natom_thiol+3, n_gold+(k-1)*natom_thiol+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 1);
    % C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 5, 2, ...
            n_gold+(k-1)*natom_thiol+3, n_gold+(k-1)*natom_thiol+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2);
    % H-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6, 4, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 1, n_gold+(k-1)*natom_thiol+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2);
  
  for i = 2:(natom_thiol_C-1)
    % C-C-C
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*(i-1)+1, 1, ...
            n_gold+(k-1)*natom_thiol+i, n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol+i+2);
    % C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*(i-1)+2, 2, ...
            n_gold+(k-1)*natom_thiol+i, n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i -1);
    % C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*(i-1)+3, 2, ...
            n_gold+(k-1)*natom_thiol+i, n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i);
    % C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*(i-1)+4, 2, ...
            n_gold+(k-1)*natom_thiol+i+2, n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i -1);
    % C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*(i-1)+5, 2, ...
            n_gold+(k-1)*natom_thiol+i+2, n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i);
    % H-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*(i-1)+6, 4, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i -1, n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i);
  end
  
  % C-C-H
  fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*natom_thiol_C-5, 3, ...
          n_gold+(k-1)*natom_thiol+natom_thiol_C, n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+ k*natom_thiol - 2);
  % C-C-H
  fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*natom_thiol_C-4, 3, ...
          n_gold+(k-1)*natom_thiol+natom_thiol_C, n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+ k*natom_thiol - 1);
  % C-C-H
  fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*natom_thiol_C-3, 3, ...
          n_gold+(k-1)*natom_thiol+natom_thiol_C, n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+ k*natom_thiol);
  % H-C-H
  fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*natom_thiol_C-2, 5, ...
          n_gold+ k*natom_thiol - 2, n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+ k*natom_thiol - 1);
  % H-C-H
  fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*natom_thiol_C-1, 5, ...
          n_gold+ k*natom_thiol - 2, n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+ k*natom_thiol);
  % H-C-H
  fprintf(fp, '     %d   %d     %d     %d     %d\n', (k-1)*n_angles_per_thiol + 6*natom_thiol_C,   5, ...
          n_gold+ k*natom_thiol - 1, n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+ k*natom_thiol);  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   PRINT DIHEDRALS  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fp, '\n');
fprintf(fp, 'Dihedrals\n');
fprintf(fp, '\n');

%%% PRINT FORMAT: DIHEDRAL_N_GOLD DIHEDRAL_TYPE ATOM_ID1 ATOM_ID2 ATOM_ID3

for k = 1: n_thiol
  
    % S-C-C-C
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 1, 4, ...
            n_gold+(k-1)*natom_thiol+1, n_gold+(k-1)*natom_thiol+2, ...
            n_gold+(k-1)*natom_thiol+3, n_gold+(k-1)*natom_thiol+4);
    % S-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 2, 5, ...
            n_gold+(k-1)*natom_thiol+1, n_gold+(k-1)*natom_thiol+2, ...
            n_gold+(k-1)*natom_thiol+3, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 3);
    % S-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 3, 5, ...
            n_gold+(k-1)*natom_thiol+1, n_gold+(k-1)*natom_thiol+2, ...
            n_gold+(k-1)*natom_thiol+3, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 4);
    % C-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 4, 2, ...
            n_gold+(k-1)*natom_thiol+4, n_gold+(k-1)*natom_thiol+3, ...
            n_gold+(k-1)*natom_thiol+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 1);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 5, 3, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 1, n_gold+(k-1)*natom_thiol+2, ...
            n_gold+(k-1)*natom_thiol+3, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 3);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 6, 3, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 1, n_gold+(k-1)*natom_thiol+2, ...
            n_gold+(k-1)*natom_thiol+3, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 4);
    % C-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 7, 2, ...
            n_gold+(k-1)*natom_thiol+4, n_gold+(k-1)*natom_thiol+3, ...
            n_gold+(k-1)*natom_thiol+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 8, 3, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2, n_gold+(k-1)*natom_thiol+2, ...
            n_gold+(k-1)*natom_thiol+3, n_gold+(k-1)*natom_thiol + natom_thiol_SC +3);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9, 3, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2, n_gold+(k-1)*natom_thiol+2, ...
            n_gold+(k-1)*natom_thiol+3, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 4);
  
  
  for i = 2:(natom_thiol_C-2)
    % C-C-C-C
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9*(i-1)+1, 1, ...
            n_gold+(k-1)*natom_thiol+i, n_gold+(k-1)*natom_thiol+i+1, ...
            n_gold+(k-1)*natom_thiol+i+2, n_gold+(k-1)*natom_thiol+i+3);
    % C-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9*(i-1)+2, 2, ...
            n_gold+(k-1)*natom_thiol+i, n_gold+(k-1)*natom_thiol+i+1, ...
            n_gold+(k-1)*natom_thiol+i+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i + 1);
    % C-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9*(i-1)+3, 2, ...
            n_gold+(k-1)*natom_thiol+i, n_gold+(k-1)*natom_thiol+i+1, ...
            n_gold+(k-1)*natom_thiol+i+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i + 2);
    % C-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9*(i-1)+4, 2, ...
            n_gold+(k-1)*natom_thiol+i+3, n_gold+(k-1)*natom_thiol+i+2, ...
            n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i -1);    
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9*(i-1)+5, 3, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i -1, n_gold+(k-1)*natom_thiol+i+1, ...
            n_gold+(k-1)*natom_thiol+i+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i + 1);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9*(i-1)+6, 3, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i -1, n_gold+(k-1)*natom_thiol+i+1, ...
            n_gold+(k-1)*natom_thiol+i+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i + 2);    
    % C-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9*(i-1)+7, 2, ...
            n_gold+(k-1)*natom_thiol+i+3, n_gold+(k-1)*natom_thiol+i+2, ...
            n_gold+(k-1)*natom_thiol+i+1, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9*(i-1)+8, 3, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i, n_gold+(k-1)*natom_thiol+i+1, ...
            n_gold+(k-1)*natom_thiol+i+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i + 1);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', (k-1)*n_dihedrals_per_thiol + 9*(i-1)+9, 3, ...
            n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i, n_gold+(k-1)*natom_thiol+i+1, ...
            n_gold+(k-1)*natom_thiol+i+2, n_gold+(k-1)*natom_thiol + natom_thiol_SC + 2*i + 2);
  end
  
    % C-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', k*n_dihedrals_per_thiol - 8, 6, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C-1, n_gold+(k-1)*natom_thiol+natom_thiol_C, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+k*natom_thiol-2);
    % C-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', k*n_dihedrals_per_thiol - 7, 6, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C-1, n_gold+(k-1)*natom_thiol+natom_thiol_C, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+k*natom_thiol-1);
    % C-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', k*n_dihedrals_per_thiol - 6, 6, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C-1, n_gold+(k-1)*natom_thiol+natom_thiol_C, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+k*natom_thiol);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', k*n_dihedrals_per_thiol - 5, 7, ...
            n_gold+k*natom_thiol-4, n_gold+(k-1)*natom_thiol+natom_thiol_C, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+k*natom_thiol-2);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', k*n_dihedrals_per_thiol - 4, 7, ...
            n_gold+k*natom_thiol-4, n_gold+(k-1)*natom_thiol+natom_thiol_C, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+k*natom_thiol-1);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', k*n_dihedrals_per_thiol - 3, 7, ...
            n_gold+k*natom_thiol-4, n_gold+(k-1)*natom_thiol+natom_thiol_C, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+k*natom_thiol);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', k*n_dihedrals_per_thiol - 2, 7, ...
            n_gold+k*natom_thiol-3, n_gold+(k-1)*natom_thiol+natom_thiol_C, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+k*natom_thiol-2);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', k*n_dihedrals_per_thiol - 1, 7, ...
            n_gold+k*natom_thiol-3, n_gold+(k-1)*natom_thiol+natom_thiol_C, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+k*natom_thiol-1);
    % H-C-C-H
    fprintf(fp, '     %d   %d     %d     %d     %d     %d\n', k*n_dihedrals_per_thiol, 7, ...
            n_gold+k*natom_thiol-3, n_gold+(k-1)*natom_thiol+natom_thiol_C, ...
            n_gold+(k-1)*natom_thiol+natom_thiol_C+1, n_gold+k*natom_thiol);
end

fclose(fp);