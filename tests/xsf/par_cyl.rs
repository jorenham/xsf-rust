#[test]
fn test_pbwa() {
    xsref::test("pbwa", "d_d-d_d", |x| xsf::pbwa(x[0], x[1]));
}

#[test]
fn test_pbdv() {
    xsref::test("pbdv", "d_d-d_d", |x| xsf::pbdv(x[0], x[1]));
}

#[test]
fn test_pbvv() {
    xsref::test("pbvv", "d_d-d_d", |x| xsf::pbvv(x[0], x[1]));
}
