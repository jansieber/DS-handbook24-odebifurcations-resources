function prob=hopf_equiv_hb(oid,funcs,xhopf,pnames,phopf,omega,n_osc,cycles)
prob=ode_isol2ep(coco_prob(),oid,funcs{:},xhopf,pnames,phopf);
[~,~,sol,data,pmat,irot] = equiv_hopf_var(prob,oid,omega,n_osc,cycles);
[prob, data] = ep_add(coco_prob(), data, sol, '-no-test', '-cache-jac');
[prob, data] = ep_add_var(prob, data, sol.var.v);
prob = ep_add_HB(prob, data, sol);
prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info('HB'));
prob=ep_HB_add_func(prob,'','sym',@v_symmetry,...
    struct('gamma',pmat,'tau',irot,'dfdxhan',funcs{2}),'inactive',...
    arrayfun(@(i){sprintf('sym%02d',i)},1:2*size(pmat,1)));
end