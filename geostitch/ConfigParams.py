# imports here


class ConfigParams:
    def __init__(
        self,
        reg_lmax: int,
        win_lmax: int,
        lon_bounds: tuple[float, float],
        lat_bounds: tuple[float, float],
        win_type: str,
        mask_mode: str = "bounds",
        mask_continents: list[str] | tuple[str, ...] | None = None,
        mask_target: str = "land",
        mask_resolution: str = "110m",
        blend_mode: str = "logistic",
        blend_lcut: int = 70,
        blend_delta: float = 8.0,
        preserve_global_low_lmax: int = 5,
        reg_noise_floor: float = 1e-11,
        glo_noise_floor: float = 1e-11,
        depth_smoothing_alpha: float = 0.0,
    ):
        """Initialize configuration parameters for merging regional and global models."""

        if reg_lmax >= 0:
            self.reg_lmax = reg_lmax
        else:
            raise ValueError("reg_lmax must be non-negative")
        if win_lmax >= 0:
            self.win_lmax = win_lmax
        else:
            raise ValueError("win_lmax must be non-negative")
        if mask_mode not in ["bounds", "continent"]:
            raise ValueError("mask_mode must be 'bounds' or 'continent'")
        self.mask_mode = mask_mode

        if self.mask_mode == "bounds":
            if len(lon_bounds) != 2:
                raise ValueError("lon_bounds should only have 2 values")
            elif lon_bounds[0] > lon_bounds[1]:
                raise ValueError("lon_bounds[0] can't be greater than lon_bounds[1]")
            else:
                for bound in lon_bounds:
                    if not (0 <= bound <= 360):
                        raise ValueError(f"lon_bounds value {bound} is not between 0 and 360")
                self.lon_min_mask, self.lon_max_mask = lon_bounds
            if len(lat_bounds) != 2:
                raise ValueError("lat_bounds should only have 2 values")
            elif lat_bounds[0] > lat_bounds[1]:
                raise ValueError("lat_bounds[0] can't be greater than lat_bounds[1]")
            else:
                for bound in lat_bounds:
                    if not (-90 <= bound <= 90):
                        raise ValueError(f"lat_bounds value {bound} is not between -90 and 90")
                self.lat_min_mask, self.lat_max_mask = lat_bounds
        else:
            self.lon_min_mask = None
            self.lon_max_mask = None
            self.lat_min_mask = None
            self.lat_max_mask = None

        if mask_target not in ["land", "ocean"]:
            raise ValueError("mask_target must be 'land' or 'ocean'")
        self.mask_target = mask_target
        self.mask_resolution = mask_resolution
        self.mask_continents = list(mask_continents) if mask_continents is not None else []

        if win_type not in ['spherical', 'rectangular']:
            raise ValueError("win_type must be 'spherical' or 'rectangular'")
        else:
            self.win_type = win_type

        # Blending controls
        if blend_mode not in ["adaptive", "logistic"]:
            raise ValueError("blend_mode must be 'adaptive' or 'logistic'")
        self.blend_mode = blend_mode

        if blend_lcut < 0:
            raise ValueError("blend_lcut must be non-negative")
        self.blend_lcut = blend_lcut

        if blend_delta <= 0:
            raise ValueError("blend_delta must be positive")
        self.blend_delta = float(blend_delta)

        if preserve_global_low_lmax < 0:
            raise ValueError("preserve_global_low_lmax must be non-negative")
        self.preserve_global_low_lmax = preserve_global_low_lmax

        if reg_noise_floor <= 0 or glo_noise_floor <= 0:
            raise ValueError("reg_noise_floor and glo_noise_floor must be positive")
        self.reg_noise_floor = float(reg_noise_floor)
        self.glo_noise_floor = float(glo_noise_floor)

        if not (0.0 <= depth_smoothing_alpha <= 1.0):
            raise ValueError("depth_smoothing_alpha must be in [0, 1]")
        self.depth_smoothing_alpha = float(depth_smoothing_alpha)
