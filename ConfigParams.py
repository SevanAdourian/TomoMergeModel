# imports here


class ConfigParams:
    def __init__(
        self,
        reg_lmax: int,
        win_lmax: int,
        lon_bounds: tuple[int],
        lat_bounds: tuple[int],
        win_type: str,
    ):
        self.reg_lmax = reg_lmax
        self.win_lmax = win_lmax
        self.lon_min_mask, self.lon_max_mask = lon_bounds
        self.lat_min_mask, self.lat_max_mask = lat_bounds
        self.win_type = win_type
