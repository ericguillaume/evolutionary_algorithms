import random
from enum import Enum
import numpy as np


GENE_ENCODING_BIT_COUNT = 4


class Gene(Enum):  # todo rename representation
    VALUE_0 = 0
    VALUE_1 = 1
    VALUE_2 = 2
    VALUE_3 = 3
    VALUE_4 = 4
    VALUE_5 = 5
    VALUE_6 = 6
    VALUE_7 = 7
    VALUE_8 = 8
    VALUE_9 = 9
    OPERATOR_PLUS = 10
    OPERATOR_MINUS = 11
    OPERATOR_TIMES = 12
    OPERATOR_DIVIDE = 13


def is_operator(gene):
    return gene.value >= Gene.OPERATOR_PLUS.value


def gene_value(gene):
    if is_operator(gene):
        raise ValueError("value: {} is not a value gene".format(gene))
    return gene.value


def apply_gene(operator, previous_value, new_value):
    if not is_operator(operator):
        raise ValueError("value: {} is not an operator gene".format(operator))
    elif operator == Gene.OPERATOR_PLUS:
        return previous_value + new_value
    elif operator == Gene.OPERATOR_MINUS:
        return previous_value - new_value
    elif operator == Gene.OPERATOR_TIMES:
        return previous_value * new_value
    elif operator == Gene.OPERATOR_DIVIDE:
        if new_value == 0:
            return np.nan
        return previous_value / new_value


def gene_to_string(gene):
    if not is_operator(gene):
        return str(gene.value)
    elif gene == Gene.OPERATOR_PLUS:
        return "+"
    elif gene == Gene.OPERATOR_MINUS:
        return "-"
    elif gene == Gene.OPERATOR_TIMES:
        return "*"
    elif gene == Gene.OPERATOR_DIVIDE:
        return "/"


class Chromosome:
    TARGET = 75.5

    @staticmethod
    def new_with_default_genes():
        binary_genes = [random.randint(0, 1) for _ in range(9 * GENE_ENCODING_BIT_COUNT)]

        return Chromosome(binary_genes)

    def __init__(self, binary_genes):
        self.binary_genes = binary_genes

    def fitness(self):
        genes = Chromosome.decode(self.binary_genes)

        if len(genes) == 0:
            return 0.0
        first_gene = genes[0]
        computed_value = first_gene.value

        current_operator = Gene.OPERATOR_PLUS
        should_gene_be_operator = True
        for idx, ge in enumerate(genes[1:]):
            # we select the first gene whose type is the one expected
            if should_gene_be_operator != is_operator(ge):
                continue

            if is_operator(ge):
                current_operator = ge
            else:
                computed_value = apply_gene(current_operator, computed_value, gene_value(ge))

            should_gene_be_operator = not should_gene_be_operator

        if computed_value == Chromosome.TARGET:
            result = 0
        elif np.isnan(computed_value):
            return float("inf")
        else:
            print("computed_value = {}".format(computed_value))
            result = 1.0 / abs(computed_value - Chromosome.TARGET + 1e6)
            print("result = {}".format(result))
        return result

    def mutate(self, mutation_rate):
        new_binary_genes = []
        for binary_gene in self.binary_genes:
            random_int = random.randint(0, int(1 / mutation_rate))
            if random_int == 0:
                print("mutation")
                new_binary_genes.append(0 if binary_gene == 1 else 1)
            else:
                new_binary_genes.append(binary_gene)
        return new_binary_genes

    def cross_over(self, other_ch):
        max_index = len(self.binary_genes)
        cross_over_idx = random.randint(0, max_index)  # todo need at least 1 gene per chromosome
        new_ch_binary_genes = self.binary_genes[0:cross_over_idx]
        new_ch_binary_genes.extend(other_ch.binary_genes[cross_over_idx:max_index])
        return Chromosome(new_ch_binary_genes)

    @staticmethod
    def encode(genes):
        binary_genes = []
        for gene in genes:
            binary_gene_seq = list(map(int, "{0:4b}".format(gene.value).replace(" ", "0")))
            assert len(binary_gene_seq) == GENE_ENCODING_BIT_COUNT
            for binary_gene in binary_gene_seq:
                assert (binary_gene is 0 or binary_gene is 1)
            binary_genes.extend(binary_gene_seq)
        return binary_genes

    @staticmethod
    def decode(binary_genes):
        genes = []
        for binary_gene_start_idx in range(0, len(binary_genes), GENE_ENCODING_BIT_COUNT):
            binary_gene = binary_genes[binary_gene_start_idx:(binary_gene_start_idx + GENE_ENCODING_BIT_COUNT)]
            if not len(binary_gene) == GENE_ENCODING_BIT_COUNT:
                print(len(binary_genes))
                print(binary_genes)
            assert (len(binary_gene) == GENE_ENCODING_BIT_COUNT)
            binary_gene_chars = list(map(str, binary_gene))
            gene_id = int("".join(binary_gene_chars), 2)
            if 0 <= gene_id <= 13:
                gene = Gene(gene_id)
                genes.append(gene)
        return genes
