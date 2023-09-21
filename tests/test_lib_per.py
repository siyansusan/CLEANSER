import io

from cleanser.lib_per import mm_counts


def test_normalize():
    mm_file = io.StringIO(
        """# 1
#2
#3
1 0 1
1 1 2
2 0 2
3 1 3"""
    )

    assert mm_counts(mm_file, 0) == {"1": 3, "2": 2, "3": 3}
    mm_file.seek(0)  # "rewind" the file
    assert mm_counts(mm_file, 1) == {"0": 3, "1": 5}
