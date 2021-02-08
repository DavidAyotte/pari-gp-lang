addhelp(torsion,{"
Ce fichier contient toutes les fonctions necessaires pour calculer les points
d'ordre fini sur une courbe elliptique E.\n
* Pour representer la courbe elliptique E : y^2 = x^3 + ax + b, on ecrit E = [a,b].
* Les points sur une courbes ellitiques correspondent a une liste de deux element P = [x,y]
* Le point a l'infinie d'une courbe elliptique est représenter par oo.

Voici la liste des fonctions dans ce fichier :\n
* polrootsrat(f) : Retourne les racines rationnelles d'un polyneme a coefficient entier
* elldbl(E, P) : retourne le point 2P
* ellptadd(E,P,Q) : retourne le point P+Q
* ellmultbyN(E, P, N) : retourne le point NP
* order(E, P) : retourne l'ordre de P
* elltorsion(E) : retourne une liste contenant tous les points de torsion rationnel de E avec leur ordre respectif.
"});

polrootsrat(f) = {
	T = 'x;
	f = Pol(f, T);
	Roots = List([]);
	A0 = polcoef(f, 0);
	while(A0 == 0,
		listput(Roots, 0);
		f = f/T;
		listsort(Roots, 1);
		if(f == 1, return(Roots));
		A0 = polcoef(f, 0);
	);
	A0 = polcoef(f, 0);
	AN = polcoef(f, poldegree(f));
	Div_AN = divisors(AN);
	Div_A0 = divisors(A0);
	for(i = 1, #Div_AN, 
		q = Div_AN[i];
		for(j = 1, #Div_A0,
			p = Div_A0[j];
			if(substpol(f, T, p/q) == 0, listput(Roots, p/q));
			if(substpol(f, T, -p/q) == 0, listput(Roots, -p/q));
		);
	);
	listsort(Roots, 1);
	return(Roots);
}
addhelp(polrootsrat, {"polrootsrat(f) \n
f : un polynome represente sous la forme d'un vecteur [a,b,c,d]. Par exemple, polynome ax^2 + bx + c s'ecrit [a,b,c].
Retourne les racines rationnelles d'un polyneme a coefficient entier f"
});

elldisc(E) = {
	return(4*(E[1])^3 + 27*(E[2])^2);
}
addhelp(elldisc, {"elldsisc(E) \n
E : une liste [a,b] representant la courbe elliptique y^2 = x^3 + ax + b.
Retourne le discriminant de E, soit la quantite : 4a^3 + 27b^2."
});

elldbl(E, P) = {
	if(P == oo, return(oo));
	a = E[1]; b = E[2];
	x = P[1]; y = P[2];
	X = (((x)^2 - a)^2 - 8*b*x)/(4*(y)^2);
	Y = ((x)^6 + 5*a*(x)^4 + 20*b*(x)^3 - 5*(a)^2*(x)^2 - 4*a*b*x - 8*(b)^2 - (a)^3)/(8*(y)^3);
	return([X, Y]);
}
addhelp(elldbl, {"elldbl(E, P)\n
E : une liste [a,b] representant la courbe elliptique y^2 = x^3 + ax + b.
P : un point sur la courbe E.
Retourne le point 2P."
});

ellptadd(E, P, Q) = {
	if(P == oo, return(Q));
	if(Q == oo, return(P));
	xP = P[1]; 	yP = P[2];
	xQ = Q[1]; 	yQ = Q[2];
	if(xP == xQ & yP == yQ, return(elldbl(E, P)));
	if(P == [xQ, -yQ], return(oo));
	L = (yQ - yP)/(xQ - xP);
	xR = (L)^2 - xP - xQ;
	yR = L*(xP - xR) - yP;
	return([xR, yR]);
}
addhelp(ellptadd, {"ellptadd(E,P,Q)\n
E : une liste [a,b] representant la courbe elliptique y^2 = x^3 + ax + b.
P, Q : deux points sur la courbes E (pas necessairement distincts).
Retourne le point P + Q."
});

ellmultbyN(E, P, N) = {
	if(N == 1, return(P));
	if(P[2] == 0, return(oo));
	NP = P;
	n = 2;
	for(n = 2, N,
		NP = ellptadd(E, P, NP);
	);
	return(NP);
}
addhelp(ellmultbyN, {"ellmultbyN(E, P, N)\n
E : une liste [a,b] representant la courbe elliptique y^2 = x^3 + ax + b.
P : un point sur la courbes E.
N : u nombre entier.
Retourne le point NP."
});

order(E, P) = {
	if(P == oo, return(1));
	for(n = 2, 10,
		if(ellmultbyN(E, P, n) == oo, return(n))
	);
	if(ellmultbyN(E, P, 12) == oo, return(12));
	return(oo);
}
addhelp(order, {"order(E, P)\n
E : une liste [a,b] representant la courbe elliptique y^2 = x^3 + ax + b.
P : un point sur la courbes E.
Retourne l'ordre de P. Si P est d'ordre infinie, la fonction retourne +oo.
"});

elltorsion(E) = {
	a = E[1]; b = E[2];
	D = elldisc(E);
	Torspts = List([[oo, 1]]);
	\\Cas y = 0 :
	f = [1, 0, a, b];
	roots = polrootsrat(f);
	for(j = 1, #roots,
		P = [roots[j], 0];
		ord = order(E, P);
		if(ord < oo, listput(Torspts, [P, ord]));
	);
	Div = divisors(D);
	for(i = 1, #Div,
		y = Div[i];
		if(D%(y^2) == 0,
			f = [1, 0, a, b-y^2];
			roots = polrootsrat(f);
			for(j = 1, #roots,
				P = [roots[j], y];
				ord = order(E, P);
				if(ord < oo, listput(Torspts, [P, ord]);listput(Torspts, [[P[1], -P[2]], ord]));
			);
		);
	);
	return(Torspts);
}
addhelp(elltorsion, {"elltorsion(E)\n
E : une liste [a,b] representant la courbe elliptique y^2 = x^3 + ax + b.
Retourne une liste contenant tous les points rationnel d'ordre fini sur la courbe E, avec leur ordre respectif.
"});

