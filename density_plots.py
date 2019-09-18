from plotnine import *
import pandas as pd


def create_density_plot(probabilities, extended_list, output_path):

    y = probabilities['gene_symbol'].map(
        lambda name: name in set(extended_list['gene_symbol']))
    x = [i for i in range(1, len(probabilities) + 1)]

    df = pd.DataFrame({'gene_symbol': x, 'Class': y})

    fig = (

        ggplot(df[df['Class'] == True], aes(x='gene_symbol', col='Class'))
        + geom_density(fill="#56B4E9", alpha = 0.3, colour = None, show_legend = True)
        + geom_rug(color='black', show_legend = True)
        + theme_bw()
        + theme(plot_title = element_text(hjust = 0.5, size=12), figure_size = (10, 6))
        + theme(legend_background = element_rect(fill="gray90", size=.5, linetype="dotted"))
        + theme(legend_title = element_text(colour="black", size=12))
        + ggtitle('Distribution of gene names from extended list')
        + xlab('Sorted in descending order gene names')
        + ylab('Density')
    )
    save_as_pdf_pages([fig])
    pass
