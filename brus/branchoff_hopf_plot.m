clear
  srun = @(i)sprintf('hopf_eq%d', i);
  hopf_eq = arrayfun(@(i){coco_bd_table(srun(i))}, 1:3);
  
  nruns = length(hopf_eq);
  names = { ...
    '$P_1$: $\Pi\sim(1\,2\,4\,3)_4$', ...
    '$P_2$: $\Pi\sim(2\,1)_2(3\,4)_1$', ...
    '$P_3$: $\Pi\sim(3\,4\,2)_1$'};
  
  proflab = {{'UST',2}, {'FP',1}, {'UST',1}};
  clri = [1,2,4];
  
  figure(1); clf;
  tl = tiledlayout(3,2);
  tl.TileSpacing = 'tight';
  nexttile([3,1]);
  ax1 = gca; hold(ax1,'on');
  
  lw  = {'linewidth', 2};
  txt = {'FontSize', 20, 'FontName', 'courier', 'FontWeight', 'bold'};
  ltx = {'Interpreter', 'latex'};
  clr = lines();
  ty  = @(r,s)strcmp(r{:,'TYPE'}, s);
  
  for i=1:nruns
    pl(i) = plot(ax1, hopf_eq{i}.lambda, hopf_eq{i}.amplitude,lw{:}, ...
      'DisplayName', names{i}, 'color', clr(clri(i),:)); %#ok<*AGROW>
    specsel{i} = ty(hopf_eq{i}, 'UST') | ty(hopf_eq{i}, 'BP') | ...
      ty(hopf_eq{i}, 'FP');
    uz = plot(ax1, hopf_eq{i}.lambda(specsel{i}), ...
      hopf_eq{i}.amplitude(specsel{i}), 'ks', lw{:});
  end
  
  hbi = find(ty(hopf_eq{2}, 'EP') & hopf_eq{2}.amplitude<1e-2);
  hb  = plot(ax1, hopf_eq{2}.lambda(hbi), hopf_eq{2}.amplitude(hbi), ...
    'kd', lw{:}, 'MarkerSize', 12, 'MarkerFaceColor', 'g', ...
    'DisplayName', 'equivariant Hopf bifurcation');
  uz.DisplayName = sprintf('UST (number of unstable$^{\\phantom{\\int}}$\nFloquet multipliers changes)');
  legend(ax1, [pl,uz,hb], ltx{:}, txt{:}, 'Location', 'NorthEast');
  xlim(ax1, [18,22]);
  ylim(ax1, [0,1]);
  xlabel(ax1, '$\epsilon\times10^{3}$', ltx{:}, txt{:})
  ylabel(ax1, 'amplitude', ltx{:}, txt{:})
  set(ax1, txt{:}, lw{:}, 'box', 'on');
  
  for i=1:nruns
    ustlabs{i} = coco_bd_labs(srun(i), proflab{i}{1});
    hopf_eqprof{i} = po_read_solution(srun(i), ustlabs{i}(proflab{i}{2}));
    nexttile(2*i); ax(i)=gca; hold(ax(i), 'on');
    t = hopf_eqprof{i}.tbp/hopf_eqprof{i}.T;
    tprofs{i} = plot(ax(i), t, hopf_eqprof{i}.xbp(:,1:2:end), lw{:});
    for k=1:4
      set(tprofs{i}(k), 'DisplayName', sprintf('$x_%d$', k));
    end
    pllab{i} = plot(ax(i), NaN, NaN, 'o', 'MarkerEdgeColor', 'k', ...
      lw{:}, 'MarkerFaceColor', clr(clri(i),:), ...
      'DisplayName', sprintf('$P_%d$', i));
    legend(ax(i), [pllab{i};tprofs{i}], txt{:}, ltx{:}, ...
      'Location','EastOutside');
    set(ax(i), txt{:}, lw{:}, 'box', 'on');
    if i<3
      ax(i).XTickLabel = {};
    end
    ax(i).YTickLabel = {};
    ax(i).YLimitMethod = 'padded';
    ax(i).XLimitMethod = 'padded';
  end
  set(tprofs{2}(4), 'LineStyle', '--');
  set(tprofs{3}(4), 'LineStyle', '--');
  set(tprofs{3}(2), 'LineStyle', ':', 'LineWidth', 4);
  xlabel(ax(3), 'time $t/T$', txt{:}, ltx{:});

%%
exportgraphics(figure(1),'../figures/equivariant_hopf.pdf')