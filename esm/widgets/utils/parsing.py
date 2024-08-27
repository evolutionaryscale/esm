def convert_range_string_to_list_of_ranges(range_str: str) -> list[tuple[int, int]]:
    def parse_range(range_str: str) -> list[tuple[int, int]]:
        result: list[tuple[int, int]] = []
        for r in range_str.split(","):
            if "-" in r:
                start, end = map(int, r.split("-"))
                result.append((start, end))
            else:
                start = end = int(r)
                result.append((start, end))
        return result

    return parse_range(range_str)
