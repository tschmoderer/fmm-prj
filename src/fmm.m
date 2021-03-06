clc
close all
clear all

% espace [0,1]x[0,1]

N = 50; % Nb de pts de discrétisations dans le sens des x et des y
dx = 1;

% matrice des temps de passage
% on padd T par des bordures (pour gérer les bords facilements) 
init = 100;
T = 10^5*ones(N+3);
T(2:end-1,2:end-1) = init*ones(N+1); % les temps effectifs de passage

% Initialisation : 
% On créé des masques AP, NB

% t = 0 les poinst acceptés sont 
ap = zeros(N+1);
ap(1,1)= 1; ap(1,2) = 1; ap(2,1) = 1; ap(2,2) = 1;
AP = zeros(N+3); 
%AP(2:end-1,2:end-1) = ap; % on applique le mm padding
%AP(9,2:end-1) = 1;
%AP(2:end-1,9) = 1;
AP(2,2:end-1) = 1; AP(16,15) = 1;
AP = logical(AP);

T(AP) = 0; % les temps de passages à l'instant 0
t = 0; % décompte du temps 

%for k = 1:10000
while sum(AP(:)) < (N+1)^2
	NB = narrow(AP);
	NB = logical(NB);

	c = velocity(N);

	% calcul du temps de passage sur la narrow band
	% On va construire TG, TD, TH, TB
	% on incruste Nb dans un truc plus gros 

	NG = zeros(size(NB)); NG(2:end-1,1:end-2) = NB(2:end-1,2:end-1); % les noeuds à gauche;
	ND = zeros(size(NB)); ND(2:end-1,3:end) = NB(2:end-1,2:end-1); % les noeuds à droite;
	NH = zeros(size(NB)); NH(1:end-2,2:end-1) = NB(2:end-1,2:end-1); % les noeuds à en haut;
	Nb = zeros(size(NB)); Nb(3:end,2:end-1) = NB(2:end-1,2:end-1); % les noeuds à bas;

	NG = logical(NG); ND = logical(ND); NH = logical(NH); Nb = logical(Nb);

	TG = T(NG); TD = T(ND); TH = T(NH); Tb = T(Nb);
	Ti = T(NB);
	ci = c(NB);

	for j = 1:length(Ti)
		if (Ti(j) < TG(j) && Ti(j) < TD(j) && Ti(j) < TH(j) && Ti(j) < Tb(j))
			error('erreur résolution de l équation');
		elseif (max(max(Ti(j) - TG(j),Ti(j) - TD(j)),0) == 0)
			if (Ti(j) - TH(j) <= Ti(j) - Tb(j))
				% (Ti-Tb)^2 = ..
				Ti(j) = Tb(j) + (dx/ci(j));
			else 
				% (Ti-TH)^2 = ..
				Ti(j) = TH(j) + (dx/ci(j));
			end
		elseif (max(max(Ti(j) - TH(j),Ti(j) - Tb(j)),0) == 0)
			if (Ti(j) - TG(j) <= Ti(j) - TD(j))
				% (Ti-TD)^2 = ..
				Ti(j) = TD(j) + (dx/ci(j));
			else 
				% (Ti-TD)^2 = ..
				Ti(j) = TG(j) + (dx/ci(j));
			end
		elseif	(Ti(j) - TG(j) <= Ti(j) - TD(j) && Ti(j) - TH(j) <= Ti(j) - Tb(j))
			% (Ti-TD)^2 + (Ti-Tb)^2=..
			delta = (2*TD(j)+2*Tb(j))^2-4*2*(TD(j)^2+Tb(j)^2-(dx/ci(j))^2);
			Ti(j) = 0.5*(TD(j) + Tb(j)) + 0.25*sqrt(delta);
		elseif (Ti(j) - TG(j) <= Ti(j) - TD(j) && Ti(j) - TH(j) >= Ti(j) - Tb(j))
			% (Ti-TD)^2 + (Ti-TH)^2=..
			delta = (2*TD(j)+2*TH(j))^2-4*2*(TD(j)^2+TH(j)^2-(dx/ci(j))^2);
			Ti(j) = 0.5*(TD(j) + TH(j)) + 0.25*sqrt(delta);
		elseif (Ti(j) - TG(j) >= Ti(j) - TD(j) && Ti(j) - TH(j) <= Ti(j) - Tb(j))
			% (Ti-TG)^2 + (Ti-Tb)^2=..
			delta = (2*TG(j)+2*Tb(j))^2-4*2*(TG(j)^2+Tb(j)^2-(dx/ci(j))^2);
			Ti(j) = 0.5*(TG(j) + Tb(j)) + 0.25*sqrt(delta);
		elseif (Ti(j) - TG(j) >= Ti(j) - TD(j) && Ti(j) - TH(j) >= Ti(j) - Tb(j))
			% (Ti-TG)^2 + (Ti-TH)^2=..
			delta = (2*TG(j)+2*TH(j))^2-4*2*(TG(j)^2+TH(j)^2-(dx/ci(j))^2);
			Ti(j) = 0.5*(TG(j) + TH(j)) + 0.25*sqrt(delta);
		end
	end

	% trouver le min de Ti 
	[m I] = min(Ti); 
    I = find(Ti == m);
	%Ti([1:length(Ti)]~=I') = init;
    Ti(~ismember(I,[1:length(Ti)])) = init;
	T(NB) = Ti;
	% mise à jour des ts AP : 
	accept = zeros(length(T(NB)),1); accept(I) = 1;
	AP(NB) = accept;
	AP = logical(AP);

	t = t+ m;

	%% function plot 

	[XX YY] = meshgrid([1:N+1],[1:N+1]);
	surf(XX,YY,T(2:end-1,2:end-1)); 
	xlabel('x')
	ylabel('y')
	zlabel('T')
	title(['time : ' num2str(t)]);
	drawnow
%	pause(0.04)
end 

