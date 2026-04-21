"""Tests for ExposoGraph.db_clients (KEGG, CTD, IARC)."""

from unittest.mock import MagicMock, patch

from ExposoGraph.db_clients.ctd import CTDClient
from ExposoGraph.db_clients.iarc import IARCClassifier, IARCGroup
from ExposoGraph.db_clients.kegg import KEGGClient, KEGGGene, KEGGPathway

# ── KEGG Client ───────────────────────────────────────────────────────────


class TestKEGGClient:
    @patch("requests.get")
    def test_get_pathway(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.text = (
            "ENTRY       hsa05204\n"
            "NAME        Chemical carcinogenesis - DNA adducts\n"
            "GENE        CYP1A1  cytochrome P450 [KO:K07408]\n"
            "            CYP1B1  cytochrome P450 [KO:K07409]\n"
            "COMPOUND    C00001\n"
        )
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = KEGGClient()
        pathway = client.get_pathway("hsa05204")

        assert isinstance(pathway, KEGGPathway)
        assert pathway.pathway_id == "hsa05204"
        assert pathway.name == "Chemical carcinogenesis - DNA adducts"
        assert "CYP1A1" in pathway.genes

    @patch("requests.get")
    def test_get_pathway_parses_numeric_gene_rows(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.text = (
            "ENTRY       hsa05204\n"
            "NAME        Chemical carcinogenesis - DNA adducts\n"
            "GENE        1543  CYP1A1; cytochrome P450 [KO:K07408]\n"
            "            1544  CYP1A2; cytochrome P450 [KO:K07409]\n"
            "            1545  CYP1B1; cytochrome P450 [KO:K07410]\n"
            "COMPOUND    C00001\n"
        )
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = KEGGClient()
        pathway = client.get_pathway("hsa05204")

        assert pathway.genes == ["CYP1A1", "CYP1A2", "CYP1B1"]

    @patch("requests.get")
    def test_get_pathway_strips_path_prefix(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.text = "ENTRY       hsa05204\nNAME        Test\n"
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = KEGGClient()
        client.get_pathway("path:hsa05204")

        url = mock_get.call_args.args[0]
        assert "path:" not in url
        assert "hsa05204" in url

    @patch("requests.get")
    def test_get_gene(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.text = (
            "ENTRY       1543\n"
            "SYMBOL      CYP1A1\n"
            "NAME        cytochrome P450 family 1 subfamily A member 1\n"
            "PATHWAY     hsa05204\n"
        )
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = KEGGClient()
        gene = client.get_gene("hsa:1543")

        assert isinstance(gene, KEGGGene)
        assert gene.symbol == "CYP1A1"
        assert "hsa05204" in gene.pathways

    @patch("requests.get")
    def test_get_gene_reads_continued_pathway_lines(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.text = (
            "ENTRY       1543\n"
            "SYMBOL      CYP1A1\n"
            "NAME        cytochrome P450 family 1 subfamily A member 1\n"
            "PATHWAY     hsa05204  Chemical carcinogenesis - DNA adducts\n"
            "            hsa00980  Metabolism of xenobiotics by cytochrome P450\n"
            "            hsa00830  Retinol metabolism\n"
        )
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = KEGGClient()
        gene = client.get_gene("hsa:1543")

        assert gene.pathways == ["hsa05204", "hsa00980", "hsa00830"]

    @patch("requests.get")
    def test_find_genes(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.text = "hsa:1543\tCYP1A1; cytochrome P450\nhsa:1545\tCYP1B1; cytochrome P450\n"
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = KEGGClient()
        results = client.find_genes("CYP1")

        assert len(results) == 2
        assert results[0]["gene_id"] == "hsa:1543"
        assert "CYP1A1" in results[0]["description"]

    @patch("requests.get")
    def test_list_pathway_genes(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.text = "path:hsa05204\thsa:1543\npath:hsa05204\thsa:1545\n"
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = KEGGClient()
        genes = client.list_pathway_genes("hsa05204")

        assert len(genes) == 2
        assert "hsa:1543" in genes

    @patch("requests.get")
    def test_empty_find_results(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.text = ""
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = KEGGClient()
        results = client.find_genes("nonexistent_gene_xyz")
        assert results == []

    def test_custom_base_url(self):
        client = KEGGClient(base_url="http://test:8080/")
        assert client.base_url == "http://test:8080"


# ── CTD Client ────────────────────────────────────────────────────────────


class TestCTDClient:
    @patch("requests.get")
    def test_get_chemical_gene_interactions(self, mock_get):
        tsv = (
            "# CTD Data\n"
            "Benzo(a)pyrene\tD001564\tCYP1A1\t1543\tHomo sapiens\tmetabolism\t\t12345678\n"
            "Benzo(a)pyrene\tD001564\tGSTM1\t2944\tHomo sapiens\tdetoxification\t\t23456789\n"
        )
        mock_resp = MagicMock()
        mock_resp.text = tsv
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = CTDClient()
        interactions = client.get_chemical_gene_interactions("Benzo(a)pyrene")

        assert len(interactions) == 2
        assert interactions[0].gene_symbol == "CYP1A1"
        assert interactions[0].chemical_name == "Benzo(a)pyrene"
        assert "12345678" in interactions[0].pubmed_ids

    @patch("requests.get")
    def test_filters_by_organism(self, mock_get):
        tsv = (
            "BaP\tD001564\tCyp1a1\t1543\tMus musculus\tmetabolism\t\t\n"
            "BaP\tD001564\tCYP1A1\t1543\tHomo sapiens\tmetabolism\t\t\n"
        )
        mock_resp = MagicMock()
        mock_resp.text = tsv
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = CTDClient()
        interactions = client.get_chemical_gene_interactions("BaP", organism="Homo sapiens")

        assert len(interactions) == 1
        assert interactions[0].organism == "Homo sapiens"

    @patch("requests.get")
    def test_empty_response(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.text = "# No results\n"
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = CTDClient()
        interactions = client.get_chemical_gene_interactions("nonexistent")
        assert interactions == []

    @patch("requests.get")
    def test_get_gene_interactions(self, mock_get):
        tsv = "BaP\tD001564\tCYP1A1\t1543\tHomo sapiens\tmetabolism\t\t\n"
        mock_resp = MagicMock()
        mock_resp.text = tsv
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        client = CTDClient()
        interactions = client.get_gene_interactions("CYP1A1")

        assert len(interactions) == 1
        assert interactions[0].gene_symbol == "CYP1A1"

    def test_custom_base_url(self):
        client = CTDClient(base_url="http://test:9090/query")
        assert client.base_url == "http://test:9090/query"


# ── IARC Classifier ───────────────────────────────────────────────────────


class TestIARCClassifier:
    def test_classify_group_1(self):
        clf = IARCClassifier()
        assert clf.classify("Benzo[a]pyrene") == IARCGroup.GROUP_1

    def test_classify_group_2a(self):
        clf = IARCClassifier()
        assert clf.classify("Acrylamide") == IARCGroup.GROUP_2A

    def test_classify_group_2b(self):
        clf = IARCClassifier()
        assert clf.classify("PhIP") == IARCGroup.GROUP_2B

    def test_classify_group_3(self):
        clf = IARCClassifier()
        assert clf.classify("Pyrene") == IARCGroup.GROUP_3

    def test_classify_unknown(self):
        clf = IARCClassifier()
        assert clf.classify("NotAChemical") is None

    def test_get_entry(self):
        clf = IARCClassifier()
        entry = clf.get_entry("BaP")
        assert entry is not None
        assert entry["group"] == "Group 1"
        assert entry["category"] == "PAH"
        assert entry["cas"] == "50-32-8"

    def test_get_entry_unknown(self):
        clf = IARCClassifier()
        assert clf.get_entry("NotAChemical") is None

    def test_list_by_group(self):
        clf = IARCClassifier()
        group_1 = clf.list_by_group(IARCGroup.GROUP_1)
        assert "Benzo[a]pyrene" in group_1
        assert "BaP" in group_1
        assert len(group_1) >= 5

    def test_list_by_category(self):
        clf = IARCClassifier()
        pahs = clf.list_by_category("PAH")
        assert "Benzo[a]pyrene" in pahs
        assert len(pahs) >= 3

    def test_all_chemicals(self):
        clf = IARCClassifier()
        assert len(clf.all_chemicals) >= 20

    def test_extra_data(self):
        clf = IARCClassifier(
            extra={"CustomChem": {"group": "Group 1", "cas": "000-00-0", "category": "Test"}},
        )
        assert clf.classify("CustomChem") == IARCGroup.GROUP_1

    def test_alias_bap(self):
        clf = IARCClassifier()
        assert clf.classify("BaP") == clf.classify("Benzo[a]pyrene")

    def test_get_entry_handles_common_case_variants(self):
        clf = IARCClassifier()
        entry = clf.get_entry("Vinyl Chloride")
        assert entry is not None
        assert entry["cas"] == "75-01-4"

    def test_short_nitrosamine_aliases_are_available(self):
        clf = IARCClassifier()
        assert clf.classify("NNK") == IARCGroup.GROUP_1
        assert clf.classify("NDMA") == IARCGroup.GROUP_2A
