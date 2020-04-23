import sys
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.lines import Line2D
import math

max_length = 401


def last_2_bits(x):
    return x & 3


def h_index_to_xy(hindex, N):
    positions = [
        [0, 0],
        [0, 1],
        [1, 1],
        [1, 0]
    ]

    tmp = positions[last_2_bits(hindex)]
    hindex = hindex >> 2

    # 2. iteratively compute coords
    x = tmp[0]
    y = tmp[1]

    n = 4
    while n <= N:

        n2 = int(n / 2)
        pos_in_small_square = last_2_bits(hindex)

        if pos_in_small_square == 0:  # lower left
            tmp = x
            x = y
            y = tmp
        elif pos_in_small_square == 1:  # upper left
            x = x
            y = y + n2
        elif pos_in_small_square == 2:  # upper right
            x = x + n2
            y = y + n2
        elif pos_in_small_square == 3:  # lower right
            tmp = y
            y = (n2 - 1) - x
            x = (n2 - 1) - tmp
            x = x + n2

        hindex = hindex >> 2
        n *= 2

    return x, y


def draw_hilbert(order, fig_width, fig_height):
    fig, ax = plt.subplots()

    # Make graph square
    fig.set_size_inches(fig_width, fig_height)

    # Move graph window a little left and down
    # scatter([-0.1],[-0.1],s=0.01)

    N = 2 ** order
    prev = (0, 0)

    print("drawing...")

    for i in range(N * N):
        curr = h_index_to_xy(i, N)
        # print(prev, curr)

        # line from prev to curr
        h_line = [prev, curr]
        (h_line_x, h_line_y) = zip(*h_line)
        ax.add_line(Line2D(h_line_x, h_line_y, linewidth=1, color='blue'))

        prev = curr

        if i % 1000 == 0:
            print(i, " done")

    plt.plot()
    plt.show()


'''
Returns the pixel co-ordinates for a hilbert curve of the given order.
Writes the pixel co-ordinates in a file
'''


def write_pixel_list_hilbert(order, file_name):
    point_list = []
    N = 2 ** order
    prev = (0, 0)

    print("-- writing --")

    pixel_count = 0
    with open(file_name, "w") as pixel_file:
        for i in range(N * N):
            curr = h_index_to_xy(i, N)
            point_list.append(curr)
            pixel_file.write(str(curr) + '\n')
            pixel_count += 1

    print("PixelCount: ", pixel_count)
    return point_list


'''
Genetically reverse an input sequence. 
A becomes T, G becomes C, and vice versa
'''


def reverse_genetically(input_seq):
    sequence = input_seq.upper()
    rev_seq = []
    for k in range(len(sequence)):
        if sequence[k] == 'A':
            rev_seq.append('T')
        elif sequence[k] == 'T':
            rev_seq.append('A')
        elif sequence[k] == 'C':
            rev_seq.append('G')
        elif sequence[k] == 'G':
            rev_seq.append('C')

    if len(rev_seq) != max_length:
        print(len(rev_seq))
    return ''.join(rev_seq)


'''
Reads the fasta file given as input. 
Removes the numbers (if any) from the sequences. 
Returns the modified list.
'''


def remove_numbers_from_sequence(file_name):
    sequence_list = []
    for record in SeqIO.parse(file_name, "fasta"):
        sequence = re.sub('[0123456789]', '', str(record.seq))
        if len(sequence) != max_length:
            print("STOP")
            break
        sequence_list.append(sequence)
    return sequence_list


'''
Converts a DNA sequence to a hilbert image array.
Writes the pixel order to the file given.
'''


def make_image(sequence, width, height, output_file_name):
    order = math.log2(width)
    point_list = write_pixel_list_hilbert(order, output_file_name)
    image_array = np.ones((width, height, 4))
    for k in range(len(sequence)):
        x, y = point_list[k]
        image_array[x][y] = mapping[sequence[k].upper()]

    return image_array


make_image("abcd", 5, 6, "hukka.fasta")