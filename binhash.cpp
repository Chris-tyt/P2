#include <string.h>
#include <iostream>

#include "zmorton.hpp"
#include "binhash.hpp"

/*@q
 * ====================================================================
 */

/*@T
 * \subsection{Spatial hashing implementation}
 *
 * In the current implementation, we assume [[HASH_DIM]] is $2^b$,
 * so that computing a bitwise of an integer with [[HASH_DIM]] extracts
 * the $b$ lowest-order bits.  We could make [[HASH_DIM]] be something
 * other than a power of two, but we would then need to compute an integer
 * modulus or something of that sort.
 *
 *@c*/

#define HASH_MASK (HASH_DIM - 1)

unsigned particle_bucket(particle_t *p, float h)
{
    unsigned ix = p->x[0] / h;
    unsigned iy = p->x[1] / h;
    unsigned iz = p->x[2] / h;
    // std::cout<<ix<<"-"<<iy<<"-"<<iz<<"xyz"<<std::endl;
    return zm_encode(ix & HASH_MASK, iy & HASH_MASK, iz & HASH_MASK);
}

unsigned particle_neighborhood(unsigned *buckets, particle_t *p, float h)
{
    /* BEGIN TASK */
    // particle_t** hash = (particle_t**) buckets;
    // int ix = p->x[0] / h;
    // int iy = p->x[1] / h;
    // int iz = p->x[2] / h;
    // // 查找相邻的所有格子
    // for (int dx = -1; dx <= 1; ++dx) {
    //     for (int dy = -1; dy <= 1; ++dy) {
    //         for (int dz = -1; dz <= 1; ++dz) {
    //             // 计算邻居格子的索引
    //             int neighbor_ix = ix + dx;
    //             int neighbor_iy = iy + dy;
    //             int neighbor_iz = iz + dz;

    //             // 确保邻居格子索引有效
    //             if (neighbor_ix < 0 || neighbor_iy < 0 || neighbor_iz < 0) continue;

    //             // 获取邻居格子的哈希索引
    //             unsigned int neighbor_hash = zm_encode(neighbor_ix& HASH_MASK, neighbor_iy& HASH_MASK, neighbor_iz& HASH_MASK);

    //             // // 遍历邻居格子中的所有粒子
    //             // particle_t* neighbor = hash[neighbor_hash];
  
    //         }
    //     }
    // }
    /* END TASK */
}

void find_neighbor_buckets(particle_t **hash, particle_t *p, float h, std::vector<unsigned>& neighbors){
    /* BEGIN TASK */
    // std::cout<<"in find"<<std::endl;

    int ix = p->x[0] / h;
    int iy = p->x[1] / h;
    int iz = p->x[2] / h;
    // std::cout<<"original:"<<ix<<"x"<<iy<<"y"<<iz<<"z"<<std::endl;

    // 查找相邻的所有格子
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                // 计算邻居格子的索引
                int neighbor_ix = ix + dx;
                int neighbor_iy = iy + dy;
                int neighbor_iz = iz + dz;
                // std::cout<<neighbor_ix<<"x"<<neighbor_iy<<"y"<<neighbor_iz<<"z"<<std::endl;

                // 确保邻居格子索引有效
                if (neighbor_ix < 0 || neighbor_iy < 0 || neighbor_iz < 0 
                  ||neighbor_ix > 20 || neighbor_iy > 20 || neighbor_iz > 20 ) continue;

                // 获取邻居格子的哈希索引
                unsigned int neighbor_hash = zm_encode(neighbor_ix& HASH_MASK, neighbor_iy& HASH_MASK, neighbor_iz& HASH_MASK);
                // std::cout<<"hash: "<<neighbor_hash<<std::endl;

                neighbors.push_back(neighbor_hash);
                // // 遍历邻居格子中的所有粒子
                // particle_t* neighbor = hash[neighbor_hash];
            }
        }
    }
    // std::cout<<"out find"<<std::endl;
    /* END TASK */
}

void hash_particles(sim_state_t *s, float h)
{
    /* BEGIN TASK */
    // std::cout<<"in hash"<<std::endl;
    // 初始化哈希表和粒子链表结构
    particle_t *p = s->part;
    particle_t **hash = s->hash;

    for (int i = 0;i<HASH_SIZE;i++){
        hash[i] = nullptr;
    }

    // 遍历每个粒子，将其插入到哈希表中
    for (int i = 0; i < s->n; ++i)
    {
        // 计算哈希索引
        unsigned int hash_index = particle_bucket(p+i, h);
        // if (hash_index >= 500){
        //     std::cout<<hash_index<<"index"<<std::endl;
        // }

        // 将粒子插入到哈希表中对应的链表中
        p[i].next = hash[hash_index];
        hash[hash_index] = &p[i];
    }
    // std::cout<<"out hash"<<std::endl;
    /* END TASK */
}
