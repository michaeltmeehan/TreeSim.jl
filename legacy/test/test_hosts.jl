
# Unit Tests
@testset "get_hosts tests" begin

    # Test 1: Single entry, no sampling, no further infections
    linelist1 = DataFrame(event = [1], child_id = [1], child_type = [1], parent_id = [0], parent_type = [0], t_birth = [1.0], t_death = [Inf], t_sam = [-1.0])
    hosts1 = get_hosts(linelist1)
    
    @test length(hosts1) == 0

    # Test 2: Single entry, with sampling
    linelist2 = DataFrame(event = [1], child_id = [1], child_type = [1], parent_id = [0], parent_type = [0], t_birth = [1.0], t_death = [Inf], t_sam = [2.0])
    hosts2 = get_hosts(linelist2)
    
    @test length(hosts2) == 2
    @test hosts2[1].leaf_times == [2.0]
    @test hosts2[1].leaf_ids == [1]
    
    # Test 3: Multiple entries, chain transmission
    linelist3 = DataFrame(event = [1, 1], child_id = [1, 2], child_type = [1, 2], parent_id = [0, 1], parent_type = [0, 1], t_birth = [1.0, 2.0], t_death = [Inf, Inf], t_sam = [2.0, -1.0])
    hosts3 = get_hosts(linelist3)
    
    @test length(hosts3) == 2
    @test hosts3[1].leaf_times == [2.0]
    @test hosts3[1].leaf_ids == [1]
    
    # Test 4: Multiple entries, with sampling and chain transmission
    linelist4 = DataFrame(event = [1, 1, 1], child_id = [1, 2, 3], child_type = [1, 2, 3], parent_id = [0, 1, 2], parent_type = [0, 1, 2], t_birth = [1.0, 2.0, 3.0], t_death = [Inf, Inf, Inf], t_sam = [4.0, -1.0, -1.0])
    hosts4 = get_hosts(linelist4)
    
    @test length(hosts4) == 2
    @test hosts4[1].leaf_times == [4.0]
    @test hosts4[1].leaf_ids == [1]
    
    # Test 5: Multiple entries, multiple roots
    linelist5 = DataFrame(event = [1, 1, 1, 1], child_id = [1, 2, 3, 4], child_type = [1, 1, 1, 1], parent_id = [0, 1, 1, 2], parent_type = [0, 1, 1, 1], t_birth = [1.0, 2.0, 3.0, 4.0], t_death = [Inf, Inf, Inf, Inf], t_sam = [5.0, -1.0, 4.5, 4.6])
    hosts5 = get_hosts(linelist5)
    
    @test length(hosts5) == 5
    @test hosts5[1].leaf_times == [2.0, 3.0, 5.0]
    @test hosts5[1].leaf_ids == [2, 3, 1]
    @test hosts5[2].leaf_times == [4.]
    @test hosts5[2].leaf_ids == [4]

end