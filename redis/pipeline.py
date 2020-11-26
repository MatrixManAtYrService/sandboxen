import conducto as co

img = co.Image(dockerfile="./Dockerfile", copy_dir=".")


def main() -> co.Serial:
    with co.Serial() as node:
        node["ping"] = co.Exec(
            "redis-cli -h redis-15233.c61.us-east-1-3.ec2.cloud.redislabs.com -p 15233 -a nO4bpNHpUne4PRearIOZrHYgU5N3wWsJ ping | grep PONG",
            image=img,
        )
    return node


if __name__ == "__main__":
    co.Image.share_directory("redis_test", ".")
    co.main(default=main)
