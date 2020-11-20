import conducto as co

IMG = co.Image()


def main() -> co.Serial:
    with co.Serial(image=IMG, requires_docker=True) as root:
        with co.Parallel(name="Init") as init:
            init["Build"] = co.Exec("sleep 3")
            init["Lint"] = co.Exec("sleep 1")
            init["Unit Test"] = co.Exec("sleep 1.5")
        root["Deploy"] = co.Exec("sleep 4")
        root["Integration Test"] = co.Exec("sleep 2")
    return root


if __name__ == "__main__":
    co.main(default=main)
