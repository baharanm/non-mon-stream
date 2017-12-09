function F = sfo_fn_sdpp(L,Y_loc,Gs_loc)
F.L_window = L;
F.Y_loc = Y_loc;
F.Gs_loc = Gs_loc;

F = class(F,'sfo_fn_sdpp',sfo_fn);
F = set(F,'current_set',-1);
