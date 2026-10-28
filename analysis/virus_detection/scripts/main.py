from .process import (
    combine_ct_dedup_and_nc2_data, 
    compute_virus_codetection_matrix, 
    write_virus_codetection_matrix, 
    write_virus_count_tables, 
    generate_matrix, 
    get_coverage_across_bam,
    build_virus_coverage_plotting_data
)
from .plot import (
    plot_virus_codetection_heatmap_discrete, plot_virus_results, plot_virus_codetection_heatmap_continuous_with_boxplot, plot_virus_codetection_with_boxplot_with_cbars,
    virus_coverage_plot_combined, virus_coverage_plot_paneled
)
from .gather_seqs import gather_seqs_from_castanet_into_fasta
from .align import make_consensus_from_input_ref
from .draw_LSV import draw_LSV_diagram, draw_LSV_coverage_plot
import os
import polars as pl

def main():
    # inputs
    cosmos_metadata = "data/COSMOS_metadata/COSMOS1-4.traps_combined_metadata.tsv"
    COSMOS_sample_config = "config/COSMOS1-4.local.json"
    nc2_data = "data/castanetv9.2-output/COSMOS1-4_castanetv9.2_nc2.tsv"
    dedup_reads_data = "data/castanetv9.2-output/COSMOS1-4_castanetv9.2_reads_dedup.tsv"
    BFV_results = "data/BFVculture-dilutions/bfv_expt_COSMOS3.csv"
    virus_acronyms = {
        "stratfordvirus": "STRV",
        "rossrivervirus": "RRV",
        "murrayvalleyencephalitisvirus": "MVEV",
        "kunjinvirus": "KUN",
        "kokoberavirus": "KOK",
        "edgehillvirus": "EHV",
        "japaneseencephalitisvirus": "JEV",
        "barmahforestvirus": "BFV",
        "sindbisvirus": "SINDV",
        "umbrevirus": "LSV",
    }
    # genome_coverage_groups = [(0, 0.05), (0.05, 0.75), (0.75, 0.95), (0.95, 1)]
    # genome_coverage_labels = ["missed", "partial", "near complete", "complete"]
    genome_coverage_groups = [(0, 0.3), (0.3, 1)]
    genome_coverage_labels = ["< 30%", '≥ 30%']
    genome_coverage_colours = [
            "#FDD49E",  # light orange
            "#FC8D59",  # orange
            "#E34A33",  # red-orange
            "#B30000",  # dark red
        ]#["#D9D9D2", "#DEDE76", "#B7C7A0", "#4E6E5D"]
    capture_positive_nc2_threshold = 0.05
    # make combined dataframe for plotting
    df = combine_ct_dedup_and_nc2_data(cosmos_metadata, nc2_data, dedup_reads_data, virus_acronyms, genome_coverage_groups, genome_coverage_labels)
    # TODO: Verify this suspected sample switch with TG. Swap only the Ct
    # profiles; nc2 and all other metadata remain assigned to their sample IDs.
    samples_to_swap = ["23-013-27955_S44", "23-013-28313_S17"]
    viruses_to_swap = ["stratfordvirus", "rossrivervirus"]
    first_sample_first_ct = df.filter((pl.col("sample") == samples_to_swap[0]) & (pl.col("virus") == viruses_to_swap[0]))["ct_value"].item()
    first_sample_second_ct = df.filter((pl.col("sample") == samples_to_swap[0]) & (pl.col("virus") == viruses_to_swap[1]))["ct_value"].item()
    second_sample_first_ct = df.filter((pl.col("sample") == samples_to_swap[1]) & (pl.col("virus") == viruses_to_swap[0]))["ct_value"].item()
    second_sample_second_ct = df.filter((pl.col("sample") == samples_to_swap[1]) & (pl.col("virus") == viruses_to_swap[1]))["ct_value"].item()
    df = df.with_columns(
        # pl.when((pl.col("sample") == samples_to_swap[0]) & (pl.col("virus") == viruses_to_swap[0]))
        # .then(second_sample_first_ct)
        pl.when((pl.col("sample") == samples_to_swap[1]) & (pl.col("virus") == viruses_to_swap[1]))
        .then(pl.lit(0))
        # .when((pl.col("sample") == samples_to_swap[0]) & (pl.col("virus") == viruses_to_swap[1]))
        # .then(second_sample_second_ct)
        .when((pl.col("sample") == samples_to_swap[1]) & (pl.col("virus") == viruses_to_swap[0]))
        .then(second_sample_second_ct)
        .otherwise(pl.col("ct_value"))
        .alias("ct_value")
    )
    # plot results
    os.makedirs("results/virus/plots", exist_ok=True)
    # # original codetection matrix
    # codetection_matrix = compute_virus_codetection_matrix(df, virus_acronyms, capture_positive_nc2_threshold)
    # categories = {
    #     "neither": {
    #         "value": 0,
    #         "colour": "#F9F7F7E4",
    #     },
    #     "PCR only": {
    #         "value": 1,
    #         "colour": "#0072B2",
    #     },
    #     "COSMOS only": {
    #         "value": 2,
    #         "colour": "#D55E00",
    #     },
    #     "PCR and COSMOS": {
    #         "value": 3,
    #         "colour": "#009E73",
    #     }
    # }
    # # plot_virus_codetection_heatmap_discrete(codetection_matrix, "results/virus/plots/virus_codetection_heatmap.png", categories)
    # # codetection coloured by ct
    # ct_matrix = generate_matrix(df, virus_acronyms, val_column="ct_value")
    # nc2_group_matrix = generate_matrix(df, virus_acronyms, val_column="genome_coverage_category")
    # outline_colour = "black"
    # nc2_boxplot_colours = {
    #     label: "#B8B8B8C1" if label == "missed" else outline_colour
    #     for label in genome_coverage_labels
    # }
    # nc2_outline_colours = {
    #     label: "None" if label == "missed" else outline_colour
    #     for label in genome_coverage_labels
    # }
    # plot_virus_codetection_heatmap_continuous_with_boxplot(
    #     df,
    #     ct_matrix,
    #     "results/virus/plots/virus_codetection_heatmap_ct.pdf", # better for rasterising of the heat map
    #     cmap="flare",
    #     outline_by=nc2_group_matrix,
    #     outline_colours=nc2_outline_colours,
    #     boxplot_colours=nc2_boxplot_colours,
    # )
    # # plot virus results
    # plot_virus_results(
    #     df,
    #     "results/virus/plots",
    #     virus_acronyms,
    #     BFV_results,
    #     genome_coverage_groups,
    #     genome_coverage_labels,
    #     genome_coverage_colours,
    # )
    # # write count and diagnostic summary tables
    # os.makedirs("results/virus/counts", exist_ok=True)
    # write_virus_codetection_matrix(
    #     df,
    #     "results/virus/counts/virus_codetection_matrix.tsv",
    #     virus_acronyms,
    #     capture_positive_nc2_threshold
    # )
    # write_virus_count_tables(
    #     df,
    #     "results/virus/counts",
    #     virus_acronyms,
    #     genome_coverage_groups,
    #     genome_coverage_labels,
    #     capture_positive_nc2_threshold,
    # )
    # # let's collect virus seqs for tree building
    viruses_to_collect = ["umbrevirus", "sindbisvirus", "murrayvalleyencephalitisvirus", "rossrivervirus", "stratfordvirus", "japaneseencephalitisvirus", "barmahforestvirus", "edgehillvirus", "kunjinvirus", "kokoberavirus"]
    # os.makedirs("results/virus/sequences", exist_ok=True)
    # for virus in viruses_to_collect:
    #     output_fasta = f"results/virus/sequences/{virus}_consensus_sequences.fasta"
    #     gather_seqs_from_castanet_into_fasta(virus, nc2_data, COSMOS_sample_config, output_fasta, virus_acronyms)
    # # make a new alignment of LSV sample to actual LSV reference for consensus calling
    # LSV_ref = "data/references/littlesussexnbi_fixed.fasta"
    # LSV_sample = "23-013-29478Cos_S17"
    # make_consensus_from_input_ref(LSV_sample, COSMOS_sample_config, LSV_ref, "results/virus/sequences", f"{LSV_sample}_consensus_on_OR077282.1")
    # # draw LSV diagram
    # alignment_file = "data/LSV/Tanyas_consensus_cov2.mafft.fasta"
    # seqs_to_draw = {
    #     "OR077282.1 |Little Sussex virus isolate AusCH19546 segment L, complete sequence":{
    #         'label': "AusCH19546 (18-02-1976)",
    #     },
    #     "23-013-29478Cos_S17_consensus_on_OR077282.1":{
    #         'label': 'COSMOS (10-03-2023)',
    #     },
    #     "MK204828.1 |Umbre virus isolate IP segment L RNA-dependent RNA polymerase gene, complete cds":{
    #         'label': "French encephalitis (28-11-2017)"
    #     },
    #     "MN395656.1 |Umbre virus segment L RNA-dependent RNA polymerase gene, complete cds":{
    #         'label': "French encephalitis (17-07-2019)"
    #     }
    # }
    # ref_seq_id = "OR077282.1 |Little Sussex virus isolate AusCH19546 segment L, complete sequence"
    # deliniation_kmers = {
    #     'segment_L': ('agtacactactaccgaagaggaaacataaac', 'aacaatttcagtcttacaaagtaaggcagca'),
    #     'segment_M': ('gagtattgaagttaaagctatataaaatagt', 'ggagttgagtccattgatattaagagatact'),
    #     'segment_S': ('agtagtgtactccctgaacctactttgtca', 'gctgttaatagtagtcagggagcacactact')
    # }
    # genes = {
    #     "nucleoprotein": {
    #         "segment": "segment_S",
    #         "start_pos": 120,
    #         "end_pos": 833,
    #         "AA_seq": "MADSKLTFEDTQDVSRSTYNPADEYDTFVGVYREFLTIDNIRIFFINAAKAKAAMAKLKTETVEVTFGTLVVKLVNNHRIGQTQREVADDELTLHRLSGYIAKYTLHLYKSAKKAGEREAITTKIVNPISAKMGFQWSVGPEIYLSTLPGTEMFLDTFKMFPLAFILLRVKRGEIKVEMAKKAMRQKYGDKDASEWLQEDIETVKGAIKTVEKLKPTLTGLAASMTAFLQQLGINTK"
    #     },
    #     "non-structural_protein": {
    #         "segment": "segment_S",
    #         "start_pos": 157,
    #         "end_pos": 396,
    #         "AA_seq": "MLVEVLITQQMNTIRLWVSTENFSQLITLESSSSTLQRPRLRWLNSRPRQLRLRLEHSSSSWLITIELDRLSVRLRTMS"
    #     },
    #     "polymerase": {
    #         "segment": "segment_L",
    #         "start_pos": 67,
    #         "end_pos": 6894,
    #         "AA_seq": "MLDQENEQALNRLAAAVARANEPAEGKEAYNAILSRRHDYFGQVFCETIGIPYRNDVRLEDILHECDQTVDTLTNKIPNITPDNFEFRDNILTIIDYKVAVSDESADYTREKYENAVDAIKHLVDFPIRVCIVRIHPHTRQVTYSDDNFMRRFAHVVYPVDFTRYTQLTEELKEKFADSPEFMDMIAYGDVTFTANWCQNGCPELYEHPIFKEFLESMPLKYQKLFEKSLQSSCYTSERWDSNLKEVKSATRDEYEDFISRSAKELYETSPDYSEPSREEILKGWEEMTERVKEQREVSSNINDQKPSAHFIWAENNPELPTANSAKICAISKLLQSINIESSFAEDFKRLGKLFDFSEDIPGYERHCEQLKRQARSNPGQVKNKKLDPKKIGNALVLWEQQFLLSAEEMGKENKLHLLKDFFGIGKHKQFAKKGLDDTETDKPTILDFNDAEVKRSAELMMQKTKLYLSGKSSLNKEHPIASEYLNEIKKCNDETYNVLNAIFQTKYWQCITDISMLMKNMLSVSQYNRNNTFRIATCANNALYGLIFPSADIKTKDASIVYCIIAVHKDEKSLFNPGALHATFKVPGGFLSISKSMRLDKKRCHRIVTSPGIFLTIATMLKGDNPTINLNDVLNFSIFASLNITKSMLSLTEPSRYMIMNSLAVSSHVREYISEKFSPYTKTLFSVHMVEKIKAGCLNAFKQREKIKTRDIYLTDIDITQKGIQDNRDLCSIWFDGNVTLKEYLNQIYIVFYLNPKALHEKHHVMIDLLKTVIEIELDQRENVPQPWSDVPKKQTVNLPILVHSLAKNLLLDTSFAHFLRSRIESRNNFRRTFATISTFTSSKSCLKTGNFHEFKIKQCQRLKKVNTSEAKKHRIANAMFVEEQDQHLEVDHANYELMRKAIPDYVDVLSTKVFDRLYELYKIGEIDDEPVIHQIFTTMRTHTKFYFTYFNKGQKTAKDREIFVGEYEAKMCMYLVERIYKERSKVNPDEMISEPGDGKLRVLERKAEQEIRFMVEKMKAKNDEINKEIEDAKDNYVLNLGKIAALQAQKYKALKLEINADMSKWSAQDVFFKYFWSVCLDPILYPEEKEHILYFFCNYMQKELVIPDEVICNILDQRKAYPNDIISQITNGLTTNHFNVKRNWLQGNFNYMSSYAHSTAMSVFKDVIKRAADLLEAEVLIHSMVHSDDNQTSVTLVDNRICPRQFTAFVINTFEKICLTFGCQANMKKTYITNHIKEFVSLFNLYGEPFSIYGRFTLTSVGDCAYIGPYEDLASRISATQTGIKHGCPASQAWLSIAIAQWMTYLTYNMLPGQQNDPQPILEIDRKQLPLELGGYISCPLSLIALAGLESGNIDFLIKLVNKYIDPRLQREPIQDQCTDIKNWDISKLSKVEIFKLKVLRYLVLDAEMDNTDIMGETSEMRGRSLLTPRKFTTAGSLRKLVSFDDFQAVQDSEQGVDSIVEYMLQNPALLVTKGETKDDYMNTVLYRYNSKKFKESLSIQNPTQLFLEQILFSHKPIVDYTGIAEKFSNLNDGTMLEDNPTIKGRMTFTDAFAMIKKDLSTLKITNEDLEIVLKFIILNDPLTVTAANAHILHVIGHSQPRLGSTACVMPEFRNLRMIKHSPALVLKCYSKNTLLLPGVQERDMERDLNHLKTFIETTKLEDKMRNRIKKNEDLKKSKDIYFEIKEVTRFYQTCYEYVKSTDYKVKVFILPAKCMTQTDFCAVLQGSLLEDKKWKLIQHLRPVTAGGHKGEVQPVVSMDFKLAREAVNLLCFFTDSFIDSRYRRDFLLKACVEMAYKDRPMVTLIEEIKKSKNRTSFLPFLYRIGELTQEDLDAYDAQKAGGRVTWNHWQSTRSLNMGKIDFQIDTHKTHLHIIGEDQQLVIAELKIPKLTANAITQAGYRLLNEKHGFRFETMNEVLVQPRKYYITYQMKSKKQYVYKIFDTATLLYERRNIKQYTGREANPIVAVCEVVPMEEEETPALSLGDIEMYNADNESVSILKTSFNESAKMRRANLSKMQYFETIGGSMLKIRFINVNALMKDTTLMSLNYDNLVNADLARFSKLLDCDGDTDSVEDGILCFSDDPMEQTTTEEIDAIPLFKVSYKISAPARLTYKNAIRELLRREITKFMDAFTFVGPSFDSPENMGFICNICCIIQLLHTNEWSTIMLNCIHICMIYLEKDIRFHKFLLNRIFFEGLNPASGKLDYKKIKEFVISLKDREMAPWNEIFTKFKQKCVTLLDKEITRIEQNTNFENILALDENTDDQSGFLLFE"
    #     },
    #     "glycoprotein precursor": {
    #         "segment": "segment_M",
    #         "start_pos": 49,
    #         "end_pos": 4446,
    #         "AA_seq": "MVIFALLLISLAGLSSQATIRCFSGGSKIMERTSNVSLSELCLKDDVSQLKIIVDHTKNDSGLFSTASVWRKWIVSDWTTCNPKKMALGTISVIDIDKDFILKSSSYICNKECTITVDKERAQVLLQSDGLNHFEVSGTTLSSGWFKLTATVTLDQTCEHVKVSCGQKSMQFHACFKNHMSCIRFLHRSIMPASFAVSICQNIELILMTTFTLIIFLLLSILTKTYICYILLPIFIPFAYAYGFIYDKTCKKCSLCGLAYHPFTKCGIHCVCGCRFENSDAMKTHREKGKCPGFKSMRTARVLCKAKGSAFTLSVLLATLALGFLTPVQGLVLDETQASYTEDQLSESFQVMQSTISIFKFATIMNAGWTFCMFLITLALYVIYKKYAHKLLNYYVIHCKECDMLHDKSGLKYNGDFTNKCGFYTCGQVEDAMGLTVHQIRNHCSANYRTRMVLNWMAVLLLATIIKDSIMLAAAAGTEFETCVKKVEMTPECTGPFMEVGICDHQQKSLKYHEVANMLLSKGTISKLDVDIIDVLPDKIAGDLQMIQSSSDFHHQQMLEYVFLLRHCDYYSSFVSDSSYSQILWQTAAKLNETNLCKHESSKKSCKCLKSNANCGESEDINTAAQTYYKNKESELGQDLTTIMYILKKMIPGSGYSYLLNATKHNKLELFKAYMTNISSKYPNNKKLKDFMAITNAILANLFALSHSSEESIPEYKQIDFSSLVGETTDMLEANFNQTISTDKITCENPKVMICISPRTFKVSSGELYVCQRSGKYYIIDTGRFKLYKLIGAGNLYCAADAQCRYGFRILTTGEMQATFKDNCKASVHTGTPGYQNDAITHCKVLASGKCKYNNQTLIDITLCANNRYYPEAKGRIFGVDPEDLCFDSACTKQQTPYNEELISACVWENPTLRPIRPRTEVSMTFQEYKDQLLKKINTDLVLHKFVRTQNLPYFIPQFKYITIQGTATTDGITSSFIFFEIPALTGASAGYKVVNNEGQELMDMVVYIKTSKTSAKYKRVYYTGPTIAINTKHQEHCIGTCPERIEHDDGWATFSKERTSNWGCEEFGCLAVETGCLFGSCQDVIRRELDVFSKLEEDKILTELCITLSHDTYCVTIDSTTPTINDYFEIQYKTVDVQTLPKMIGVTHHKLLRGQINELGSFGKFCGNIQVVGNQTYGQAAVKFDYICHAAQRKDVVIRKCLENNFQSCKLLTHDTDLILEESETSVEVINNKRITGQIALKVLLGDFNYKQYTANVEADIEANCVGCTNCIKGITCEIKVKTEIEASCKVEPPCPSYTNRIIIKPGRDTHTLLMKCDKKITSGELTIKICNMQIQAHLSTIESNDQLELSTGDQSTYIHEEDLRCSTWICKVRDEGISFIFKPLTDWLGSFTKPILIAGCIIIVAIILIYILMPMCARLRDLLKKNEYEYLQDTKVIVKPKATASSIGLRSKKINIQSSKAVD"
    #     }
    # }
    # draw_LSV_diagram(alignment_file, seqs_to_draw, ref_seq_id, deliniation_kmers, genes, "results/virus/sequences/LSV_diagram.pdf")
    # # draw LSV coverage plot
    # LSV_bam = os.path.join("results/virus/sequences", f"{LSV_sample}_consensus_on_OR077282.1.bam")
    # LSV_coverage_df = get_coverage_across_bam(LSV_bam, "results/virus/sequences/LSV_coverage.tsv")
    # draw_LSV_coverage_plot(LSV_coverage_df, "results/virus/sequences/LSV_coverage_plot.pdf")

    # updated and hopefully final version of codetectin matrix with cbars and coloured by nc2
    # metadata_df = pl.read_csv(cosmos_metadata, separator="\t")
    # COASTAL_THRESHOLD_KM = 20
    # URBAN_CODES = [1, 2, 3, 4, 5]
    # FRINGE_CODES = [6, 7]
    # RURAL_CODES = [8, 9, 10, 11]
    # URBAN_RURAL_HEATMAP_ORDER = ["Rural", "Met_Fringe", "Urban"]
    # # if >20km from coast, is coastal
    # metadata_df = metadata_df.with_columns(
    #     pl.when(pl.col("distance_to_coast_km") > COASTAL_THRESHOLD_KM)
    #     .then(pl.lit(0))
    #     .otherwise(pl.lit(1))
    #     .alias("coastal"),
    #     pl.when(pl.col("council_olg_code").is_in(URBAN_CODES))
    #     .then(pl.lit("Urban"))
    #     .when(pl.col("council_olg_code").is_in(FRINGE_CODES))
    #     .then(pl.lit("Met_Fringe"))
    #     .when(pl.col("council_olg_code").is_in(RURAL_CODES))
    #     .then(pl.lit("Rural"))
    #     .otherwise(None)
    #     .alias("urban_rural"),
    # )
    # plot_virus_codetection_with_boxplot_with_cbars(
    #     df,
    #     metadata_df,
    #     plot_feature="nc2",
    #     cbar_features=["coastal", "urban_rural"],
    #     boxplot_x_feature="genome_coverage_category",
    #     boxplot_y_feature="ct_value",
    #     boxplot_order=genome_coverage_labels,
    #     virus_acronyms=virus_acronyms,
    #     cmap="viridis_r",
    #     output_path="results/virus/plots/virus_codetection_heatmap_ct_with_cbars.png",
    # )

    # let's get coverage across the viruses
    df = build_virus_coverage_plotting_data(viruses_to_collect, nc2_data, COSMOS_sample_config, virus_acronyms, detection_threshold=0.05, metadata_file=cosmos_metadata) # if change must re run
    df.write_csv("results/virus/plots/virus_coverage_plotting_data.csv")
    df = pl.read_csv("results/virus/plots/virus_coverage_plotting_data.csv")
    virus_coverage_plot_output = "results/virus/plots/virus_coverage_plot_combined.png"
    virus_coverage_plot_combined(df, virus_coverage_plot_output)
    virus_coverage_plot_output_paneled = "results/virus/plots/virus_coverage_plot_paneled.png"
    virus_coverage_plot_paneled(df, virus_coverage_plot_output_paneled)

if __name__ == "__main__":
    main()
