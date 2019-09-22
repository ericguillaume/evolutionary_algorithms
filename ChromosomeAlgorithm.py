from Chromosome import Chromosome, Gene, gene_to_string


# solution: 5 / 2 + 9 * 7 - 5


class ChromosomeAlgorithm:
    NB_CHROMOSOMES = 5
    MUTATION_RATE = 1e-3

    def __init__(self):
        self.chromosomes = [Chromosome.new_with_default_genes() for _ in range(ChromosomeAlgorithm.NB_CHROMOSOMES)]

    def step(self):
        new_chromosomes = []
        ch_with_fitness = [[ch.fitness(), ch] for ch in self.chromosomes]
        ch_with_fitness = sorted(ch_with_fitness, key=lambda x: x[0])

        assert (ChromosomeAlgorithm.NB_CHROMOSOMES >= 2)
        best_ch = ch_with_fitness[0][1]
        second_best_ch = ch_with_fitness[1][1]

        best_ch_fitness = ch_with_fitness[0][0]
        print("best_ch_fitness = {}".format(best_ch_fitness))
        if best_ch_fitness == 0:
            return best_ch

        new_chromosomes = list(map(lambda x: x[1], ch_with_fitness[0:int(ChromosomeAlgorithm.NB_CHROMOSOMES / 2)]))

        for _ in range(int(ChromosomeAlgorithm.NB_CHROMOSOMES / 2), ChromosomeAlgorithm.NB_CHROMOSOMES):
            new_ch = best_ch.cross_over(second_best_ch)
            new_chromosomes.append(new_ch)

        # mutate
        for ch in new_chromosomes:
            ch.mutate(ChromosomeAlgorithm.MUTATION_RATE)

        assert (len(new_chromosomes) == ChromosomeAlgorithm.NB_CHROMOSOMES)

        self.chromosomes = new_chromosomes

        return best_ch


def run():
    # encoded = Chromosome.encode([Gene.VALUE_0,
    #                    Gene.VALUE_0,
    #                    Gene.VALUE_0,
    #                    Gene.VALUE_0,
    #                    Gene.VALUE_0,
    #                    Gene.VALUE_0,
    #                    Gene.VALUE_0,
    #                    Gene.VALUE_0,
    #                    Gene.VALUE_0])
    # decoded = Chromosome.decode(encoded)
    # print(decoded)


    algo = ChromosomeAlgorithm()
    best_ch = None
    while best_ch is None:
        print("next iteration")
        best_ch = algo.step()
        break

    decoded = Chromosome.decode(best_ch.binary_genes)
    solution_gene_str = map(gene_to_string, decoded)
    solution_str = " ".join(solution_gene_str)  # pour print le ch mettre dans ch
    print("best_ch = {}".format(solution_str))


if __name__ == "__main__":
    run()
