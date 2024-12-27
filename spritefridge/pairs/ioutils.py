import pysam as ps
import itertools as it

from collections import defaultdict, namedtuple

Position = namedtuple('Position', ['chrom', 'start', 'end', 'strand'])
def read_to_pos(read, sep):
    bc_string = read.query_name.split(sep)[-1]
    pos = Position(
        chrom = read.reference_name, 
        start = read.reference_start,
        end = read.reference_end,
        strand = '+' if not read.is_reverse else '-'
    )
    return bc_string, pos


def read_bam(bamfile, bcsep):
    with ps.AlignmentFile(bamfile, 'rb') as bam:
        # until_eof relies on filtering unmapped reads beforehand
        for read in bam.fetch(until_eof = True):
            yield read_to_pos(read, bcsep)


def clusters_by_size(clusters):
    sized_clusters = defaultdict(list)
    for k, cluster in clusters.items():
        sized_clusters[len(cluster)].append([k, cluster])

    return sized_clusters


def write_pairs(clusters_of_size, filename, size):
    with open(filename, 'w') as file:
        for i, (k, cluster) in enumerate(clusters_of_size):
            for positions in it.combinations(cluster, 2):
                # ensures upper triangle
                pos1, pos2 = sorted(positions)
                #columns: readID chr1 pos1 chr2 pos2 strand1 strand2
                pair = [
                    # +1 because pairs is 1-based
                    pos1.chrom, pos1.start + 1, 
                    pos2.chrom, pos2.start + 1, 
                    pos1.strand, pos2.strand
                ]
                file.write(
                    '\t'.join([f'c_{size}_{i}', *[str(p) for p in pair]]) + '\n'
                )


def clusters_to_bed(size, clusters_of_size, file):
    for i, (k, cluster) in enumerate(clusters_of_size):
        for pos in cluster:
            file.write(
                '\t'.join(
                    [str(pos.chrom), str(pos.start), str(pos.end), f'c_{size}_{i}']
                ) + '\n'
            )


def write_bed(sized_clusters, filename, min_c_size, max_c_size):
    with open(filename, 'w') as file:
        for size, clusters_of_size in sized_clusters.items():
            if size < min_c_size or size > max_c_size:
                continue

            clusters_to_bed(
                size,
                clusters_of_size,
                file
            )


def clusters_to_pairs(clusters, fileprefix, min_c_size, max_c_size):
    sized_clusters = clusters_by_size(clusters)
    stats = {}
    for size, clusters_of_size in sized_clusters.items():
        if size < min_c_size or size > max_c_size:
            continue
        
        stats[size] = len(clusters_of_size)
        write_pairs(
            clusters_of_size, 
            f'{fileprefix}_{size}.pairs',
            size
        )
    
    write_bed(
        sized_clusters, 
        f'{fileprefix}.bed',
        min_c_size,
        max_c_size
    )

    return stats


def write_stats(stats, filename):
    with open(filename, 'w') as file:
        file.write('clustersize\tnumber_of_clusters\n')
        for size, num in stats.items():
            file.write(
                f'{size}\t{num}\n'
            )
            