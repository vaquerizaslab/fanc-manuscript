import os
import fanc
import fanc.plotting as fancplot
import fanc.plotting.statistics as fanc_stats_plot
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from fanc.plotting.base_plotter import GenomeCoordFormatter
import pickle

output_folder = './'

genome_file = 'hg19.ra.fa'
tads_file = 'rao2014_gm12878_arrowhead_domains_hg19.bed'
loops_file = 'rao2014_gm12878_hiccups_loops_hg19.bedpe'
insulation_file = 'rao2014_gm12878_10kb.ii'
directionality_file = 'rao2014_gm12878_10kb.di'
genes_file = 'gencode.v19.annotation_sorted.gtf.gz'
ctcf_file = 'stanford_ctcf_nc_FE.bw'
ctcf = fanc.load(ctcf_file)

hic_10kb_file = 'rao2014_gm12878_10kb.hic'
hic_10kb = fanc.load(hic_10kb_file)

hic_500kb_file = 'rao2014_gm12878_500kb.hic'
hic_500kb = fanc.load(hic_500kb_file)

ab_output_file = os.path.join(output_folder, 'hic_500kb.ab')
ab_enrichment_output_file = os.path.join(output_folder, 'hic_500kb.ab_enrichment.npy')
tad_am_output_file = os.path.join(output_folder, 'hic_10kb.tads.agg')
loop_am_output_file = os.path.join(output_folder, 'hic_10kb.loops.agg')


plotting_region = fanc.GenomicRegion.from_string('chr1:6mb-7mb')
triangular_plotting_region = fanc.GenomicRegion.from_string('chr18:12.5mb-13.5mb')

try:
    os.mkdir(output_folder)
except OSError:
    pass

# set up figure
matplotlib.rcParams.update({'font.size': 7})

fig = plt.figure(figsize=(8, 8), dpi=300)
gs = GridSpec(120, 120, hspace=0.2, wspace=0.2)

# Set up axes
# left plotting area
left_large_row_height, left_small_row_height, left_panel_row_gutter, left_panel_row_gutter_small = 16, 8, 8, 4
left_large_col_width, left_small_col_width, left_panel_col_gutter = 16, 16, 16
left_panel_col_colorbar, left_panel_col_colorbar_gutter_left, left_panel_col_colorbar_gutter_right = 2, 2, 12

row1 = slice(0, left_large_row_height, 1)
row2 = slice(row1.stop + left_panel_row_gutter, row1.stop + left_panel_row_gutter + left_large_row_height, 1)
row3 = slice(row2.stop + left_panel_row_gutter, row2.stop + left_panel_row_gutter + left_large_row_height, 1)
row4 = slice(row3.stop + left_panel_row_gutter_small, row3.stop + left_panel_row_gutter_small + left_small_row_height, 1)
row5 = slice(row4.stop + left_panel_row_gutter, row4.stop + left_panel_row_gutter + left_large_row_height, 1)

col1 = slice(0, left_small_col_width, 1)
col2 = slice(col1.stop + left_panel_col_colorbar_gutter_left, col1.stop + left_panel_col_colorbar_gutter_left + left_panel_col_colorbar, 1)
col3 = slice(col2.stop + left_panel_col_colorbar_gutter_right, col2.stop + left_panel_col_colorbar_gutter_right + left_large_col_width, 1)
col4 = slice(col3.stop + left_panel_col_colorbar_gutter_left, col3.stop + left_panel_col_colorbar_gutter_left + left_panel_col_colorbar, 1)

ax_hic = plt.subplot(gs[row1, col3])
cax_hic = plt.subplot(gs[row1, col4])
ax_expected = plt.subplot(gs[row2, col1])
ax_oe = plt.subplot(gs[row2, col3])
cax_oe = plt.subplot(gs[row2, col4])
ax_ab = plt.subplot(gs[row3, col1])
cax_ab = plt.subplot(gs[row3, col2])
ev_ax = plt.subplot(gs[row4, col1])
ax_ab_enrichment = plt.subplot(gs[row3, col3])
cax_ab_enrichment = plt.subplot(gs[row3, col4])
barplot_ax = plt.subplot(gs[row4, col3])
ax_aggregate_tads = plt.subplot(gs[row5, col1])
ax_aggregate_loops = plt.subplot(gs[row5, col3])
cax_aggregate_loops = plt.subplot(gs[row5, col4])

# right plotting area
right_column_start, right_column_width = col4.stop + left_panel_col_gutter + 15, 30
right_column_large_height, right_column_small_height = 10, 6
row_right1 = slice(0, right_column_large_height, 1)
row_right2 = slice(row_right1.stop + left_panel_row_gutter, row_right1.stop + left_panel_row_gutter + right_column_small_height, 1)
row_right3 = slice(row_right2.stop + left_panel_row_gutter_small, row_right2.stop + left_panel_row_gutter_small + right_column_large_height, 1)
row_right4 = slice(row_right3.stop + left_panel_row_gutter_small, row_right3.stop + left_panel_row_gutter_small + right_column_small_height, 1)
row_right5 = slice(row_right4.stop + left_panel_row_gutter, row_right4.stop + left_panel_row_gutter + right_column_large_height, 1)
row_right6 = slice(row_right5.stop + left_panel_row_gutter_small, row_right5.stop + left_panel_row_gutter_small + right_column_small_height, 1)
row_right7 = slice(row_right6.stop + left_panel_row_gutter, row_right6.stop + left_panel_row_gutter + right_column_large_height, 1)

gene_row = slice(row_right6.stop + left_panel_row_gutter, row5.stop)

col_right1 = slice(right_column_start, right_column_start + right_column_width, 1)
col_right2 = slice(right_column_start + right_column_width + 2, right_column_start + right_column_width + 4, 1)

ax_hic_triangular = plt.subplot(gs[row_right1, col_right1])
cax_hic_triangular = plt.subplot(gs[row_right1, col_right2])
ax_ctcf = plt.subplot(gs[row_right2, col_right1])
ax_ins_array = plt.subplot(gs[row_right3, col_right1])
cax_ins_array = plt.subplot(gs[row_right3, col_right2])
ax_ins = plt.subplot(gs[row_right4, col_right1])
ax_dir_array = plt.subplot(gs[row_right5, col_right1])
cax_dir_array = plt.subplot(gs[row_right5, col_right2])
ax_dir = plt.subplot(gs[row_right6, col_right1])
ax_genes = plt.subplot(gs[gene_row, col_right1])


# 1. Hi-C square plot
p_hic = fancplot.SquareMatrixPlot(hic_10kb, norm='lin', vmin=0.0, vmax=0.03, ax=ax_hic, cax=cax_hic,
                                  draw_tick_legend=False, draw_minor_ticks=False)
p_hic.plot(plotting_region)

# due to limited space we only want labels at the start and end of the plotting window
ax_hic.set_xticks([plotting_region.start, plotting_region.end])
ax_hic.set_yticks([plotting_region.start, plotting_region.end])

# similarly, we only want to label the min and max values on the colorbar
p_hic.colorbar.set_ticks([0, 0.03])
p_hic.colorbar.set_label("Normalised contact\nprobability")


# 2. Distance decay plot
fanc_stats_plot.distance_decay_plot(hic_10kb, chromosome='chr1', ax=ax_expected, tight=False)


# 3. Hi-C O/E plot
p_oe = fancplot.SquareMatrixPlot(hic_10kb, oe=True, log=True, norm='lin', vmin=-2, vmax=2, ax=ax_oe, cax=cax_oe,
                                 colormap='RdBu_r', draw_tick_legend=False, draw_minor_ticks=False)
p_oe.plot(plotting_region)

# adjust label locations for neater plot
ax_oe.set_xticks([plotting_region.start, plotting_region.end])
ax_oe.set_yticks([plotting_region.start, plotting_region.end])
p_oe.colorbar.set_ticks([-2, 0, 2])
p_oe.colorbar.set_label("Log2-O/E\ncontacts")

# 4. AB compartment plot
if not os.path.exists(ab_output_file):
    ab = fanc.ABCompartmentMatrix.from_hic(hic_500kb, file_name=ab_output_file)
else:
    ab = fanc.ABCompartmentMatrix(ab_output_file, mode='r')
p_ab = fancplot.SquareMatrixPlot(ab, norm='lin', vmin=-1, vmax=1, colormap='bwr', ax=ax_ab, cax=cax_ab,
                                 draw_tick_legend=False, draw_minor_ticks=False)
p_ab.plot('chr1')

# adjust label locations for neater plot
ax_ab.set_xticks([1, 249000000])
ax_ab.set_xticklabels([])
ax_ab.set_yticks([1, 249000000])
p_ab.colorbar.set_ticks([-1, 0, 1])
p_ab.colorbar.set_label("Pearson correlation")

# 5. AB saddle plot (500kb?)
if os.path.exists(ab_enrichment_output_file):
    with open(ab_enrichment_output_file, 'rb') as f:
        ab_enrichment_matrix, ev, cutoffs = pickle.load(f)
else:
    ev = ab.eigenvector(genome=genome_file, exclude_chromosomes=['chrY'])
    ab_enrichment_matrix, cutoffs = ab.enrichment_profile(hic_500kb,
                                                          percentiles=list(range(0, 100, 5)),
                                                          eigenvector=ev)
    with open(ab_enrichment_output_file, 'wb') as o:
        pickle.dump([ab_enrichment_matrix, ev, cutoffs], o)

distances = [r.start for r in ab.regions('chr1')]

ev_ax.axhline(0, color='#cccccc', linestyle='--')
ev_ax.plot(distances, ev[:len(distances)], linewidth=0.5, color='grey')

ev_ax.set_xlim(ax_ab.get_xlim())
# Demo how to use the axis formatter to display genome coordinates
ev_ax.xaxis.set_major_formatter(GenomeCoordFormatter('chr1',
                                                     minor_div=5,
                                                     display_chromosome=True,
                                                     display_scale=False))
ev_ax.spines['right'].set_visible(False)
ev_ax.spines['top'].set_visible(False)
ev_ax.set_xticks([1, 249000000])
ev_ax.set_ylabel('EV')

fanc_stats_plot.saddle_plot(ab_enrichment_matrix, cutoffs, colormap='RdBu_r',
                            vmin=-2, vmax=2, only_gc=False, fig=fig,
                            axes=[ax_ab_enrichment, barplot_ax, cax_ab_enrichment])
# nicer tick placement than default:
barplot_ax.set_yticks([-0.2, 0, 0.2])


# 6. TAD aggregate plot
tads = fanc.load(tads_file)
if not os.path.exists(tad_am_output_file):
    tad_am = fanc.AggregateMatrix.from_regions(hic_10kb, tads.regions, log=True, oe=True,
                                               file_name=tad_am_output_file, cache=False)
else:
    tad_am = fanc.AggregateMatrix(tad_am_output_file)

fanc_stats_plot.aggregate_plot(tad_am, labels=None, vmin=-1, vmax=1,
                               oe=False, log=False, colormap='bwr', ax=ax_aggregate_tads,
                               relative_label_locations=(1/3, 2/3), plot_colorbar=False)


# 7. Loop aggregate plot
loops = fanc.load(loops_file)
if not os.path.exists(loop_am_output_file):
    loop_am = fanc.AggregateMatrix.from_center_pairs(hic_10kb, loops, log=True, oe=True,
                                                     file_name=loop_am_output_file, cache=False)
else:
    loop_am = fanc.AggregateMatrix(loop_am_output_file)

fanc_stats_plot.aggregate_plot(loop_am, labels=('loop anchor',), vmin=-1, vmax=1,
                               oe=False, log=False, colormap='bwr', ax=ax_aggregate_loops,
                               relative_label_locations=(0.5,), cax=cax_aggregate_loops)
ax_aggregate_loops.set_yticklabels([])

# 8. Hi-C plot triangular
p_hic_triangular = fancplot.HicPlot(hic_10kb, norm='lin', ax=ax_hic_triangular, cax=cax_hic_triangular,
                                    max_dist=750000,
                                    vmin=0, vmax=0.03, draw_tick_legend=False, draw_minor_ticks=False)
p_hic_triangular.plot(triangular_plotting_region)

ax_hic_triangular.set_xticks([triangular_plotting_region.start,
                              triangular_plotting_region.center,
                              triangular_plotting_region.end])
p_hic_triangular.colorbar.set_ticks([0, 0.03])
p_hic_triangular.colorbar.set_label("Normalised\ncontact probability")

# 9. Insulation array
insulation = fanc.load(insulation_file)

p_ins_array = fancplot.GenomicVectorArrayPlot(insulation, y_scale='log', ax=ax_ins_array,
                                              colormap='bwr', cax=cax_ins_array,
                                              vmin=-1, vmax=1, draw_tick_legend=False,
                                              draw_minor_ticks=False, draw_tick_labels=False)
p_ins_array.plot(triangular_plotting_region)
p_ins_array.colorbar.set_label('Insulation score')

ax_ins_array.set_xticks([triangular_plotting_region.start,
                         triangular_plotting_region.center,
                         triangular_plotting_region.end])
ax_ins_array.set_yticks([2, 10, 100])
ax_ins_array.set_yticklabels(['20kb', '100kb', '1mb'])
ax_ins_array.set_ylabel('Window\nsize (bp)')

# 10. Insulation score
insulation_single = insulation.score_regions(10)

p_ins = fancplot.LinePlot(insulation_single, ax=ax_ins, style='mid', plot_kwargs={'color': '#79C7C5'},
                          draw_tick_legend=False, draw_minor_ticks=False)
p_ins.plot(triangular_plotting_region)
ax_ins.set_xticks([triangular_plotting_region.start,
                   triangular_plotting_region.center,
                   triangular_plotting_region.end])
ax_ins.set_ylabel('Insulation\nscore 100kb')

# 11. Directionality flame
directionality = fanc.load(directionality_file)


p_dir_array = fancplot.GenomicVectorArrayPlot(directionality, y_scale='log', ax=ax_dir_array,
                                              colormap='PuOr', cax=cax_dir_array,
                                              vmin=-0.1, vmax=0.1,
                                              draw_tick_legend=False,
                                              draw_minor_ticks=False, draw_tick_labels=False)
p_dir_array.plot(triangular_plotting_region)
p_dir_array.colorbar.set_label('Directionality index')

ax_dir_array.set_xticks([triangular_plotting_region.start,
                         triangular_plotting_region.center,
                         triangular_plotting_region.end])
ax_dir_array.set_ylabel('Window\nsize (bp)')

# 12. Directionality score
directionality_single = directionality.score_regions(1000000)

p_dir = fancplot.LinePlot(directionality_single, ax=ax_dir, style='mid', plot_kwargs={'color': '#79C7C5'},
                          draw_tick_legend=False, draw_minor_ticks=False)
p_dir.plot(triangular_plotting_region)
ax_dir.set_xticks([triangular_plotting_region.start,
                   triangular_plotting_region.center,
                   triangular_plotting_region.end])
ax_dir.set_ylabel('Directionality\nindex 1Mb')

# 13. Genes
genes = fanc.load(genes_file)
p_genes = fancplot.GenePlot(genes, feature_types=('exon',),
                            color_score=False, box_height=0.9, show_labels=True,
                            label_field='gene_name', font_size=7, arrow_size=0, show_arrows=False,
                            line_width=1, group_by='gene_id', draw_tick_legend=False,
                            draw_tick_labels=True, draw_minor_ticks=False,
                            ax=ax_genes, squash=True, min_gene_size=3000,
                            collapse=False,
                            include={
                                'gene_type': {'processed_transcript', 'rRNA',
                                              'protein_coding', 'IG_V_gene'}
                            })
p_genes.plot(triangular_plotting_region)
ax_genes.set_xticks([triangular_plotting_region.start,
                     triangular_plotting_region.center,
                     triangular_plotting_region.end])

# 14. CTCF enrichment
p_ctcf = fancplot.LinePlot(ctcf, ax=ax_ctcf, style='mid', plot_kwargs={'color': 'black'}, fill=False,
                           draw_tick_labels=False)
p_ctcf.plot(triangular_plotting_region)
ax_ctcf.set_xticks([triangular_plotting_region.start,
                   triangular_plotting_region.center,
                   triangular_plotting_region.end])
ax_ctcf.set_ylabel('CTCF occupancy\n(fold-enrichment\nover input)')

fig.savefig(os.path.join(output_folder, "figure3.png"))
plt.close(fig)
