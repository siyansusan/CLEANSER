import io

from cleanser.normalize_mRNA_lib_size import normalize


def test_normalize():
    cell_values = io.StringIO(
        """a\t1
b\t2
c\t3
d\t4"""
    )

    assert normalize(cell_values) == [("a", 0.5), ("b", 1.0), ("c", 1.5), ("d", 2.0)]
