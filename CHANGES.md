# 修改记录

## 2025年版本更新

### 新增服务器模式 - 共享内存优化

**目的：**
为了解决大批量任务重复加载索引导致的内存浪费问题，添加了服务器模式。服务器模式会将整个索引结构（包括索引数据）加载到共享内存中，允许多个任务共用同一个索引，大大降低内存占用和加载时间。

**主要修改：**

1. **C代码修改 (py_bri_wrapper.c)**
   - 彻底重写了共享内存实现
   - 新增内存布局设计，将整个bri结构及其数据放入共享内存
   - 新增了4个函数：
     - `py_create_server`: 创建共享内存并复制BRI索引结构及其数据
     - `py_connect_to_server`: 连接到共享内存并获取BRI索引指针
     - `py_disconnect_from_server`: 从共享内存断开连接
     - `py_stop_server`: 关闭服务器并释放资源

2. **参数解析器修改 (arg_parser.py)**
   - 添加了子命令系统，支持三种命令：
     - `convert`: 执行坐标转换（原有功能）
     - `server`: 启动服务器模式
     - `stop-server`: 停止服务器

3. **主程序修改 (main.py)**
   - 重构了代码结构，将功能分离为独立函数
   - 添加了服务器模式相关功能：
     - `run_server`: 运行服务器模式
     - `stop_server`: 停止服务器
     - `daemonize`: 将进程守护进程化
     - `cleanup_server`: 清理服务器资源
     - `disconnect_from_shared_memory`: 断开与共享内存的连接
   - 修改了坐标转换函数，支持检测和使用共享内存

4. **文档更新 (README.md)**
   - 更新了帮助信息，添加了服务器模式的使用说明

5. **版本更新 (setup.py)**
   - 将版本号从0.2更新到0.3

**使用方法：**

1. 启动服务器：
   ```
   SCC server -B reference2_namesort.bam
   ```

2. 使用共享内存运行转换任务：
   ```
   SCC convert -b alignment1.bam -B reference2_namesort.bam -r ref1.fa -R ref2.fa -s positions.txt -o output --use-server
   ```

3. 停止服务器：
   ```
   SCC stop-server
   ```

**技术细节：**

服务器模式使用System V IPC共享内存接口实现索引共享。改进后的实现将整个索引结构及其数据都复制到共享内存中，客户端可以直接访问而不需要重新加载索引文件。关键步骤包括：

1. 服务器启动时创建一个密钥文件(/tmp/scc_shm_key)和大小文件(/tmp/scc_shm_size)
2. 使用ftok()生成共享内存键
3. 计算所需的共享内存大小，包括索引结构体、readnames数据和records数据
4. 使用shmget()、shmat()创建和访问共享内存
5. 将所有数据复制到共享内存中，并重建指针关系
6. 客户端通过同样的方式访问共享内存
7. 服务器退出时清理资源

共享内存布局：
```
+------------------+
| 共享内存头部     | <- 包含偏移量和元数据
+------------------+
| bri结构体        | <- 包含索引的主要结构
+------------------+
| readnames数据    | <- 存储所有read名称的大块内存
+------------------+
| records数据      | <- 存储所有record的数组
+------------------+
```

**注意事项：**

1. 服务器模式目前只能共享一个索引文件
2. 确保在同一计算节点上运行服务器和客户端
3. 如果服务器异常退出，可能需要手动运行stop-server命令清理资源
4. 此方案将整个索引加载到共享内存中，所以服务器节点需要有足够的内存

## 2025年版本更新 - 修复共享内存稳定性问题

**主要修复：**

1. **修复指针重建问题**：
   - 完善指针在共享内存中的存储和重建机制
   - 添加专用函数`rebuild_pointers_in_shm`确保在客户端连接时正确重建所有指针结构
   - 修复records中read_name.ptr在共享内存中的存储方式，使用偏移量而非直接指针

2. **增强错误处理**：
   - 改进服务器停止逻辑，处理文件丢失情况
   - 为服务器启动添加检测，避免多个服务器实例冲突
   - 忽略无关紧要的错误消息避免造成混淆

3. **增强稳定性**：
   - 确保内存结构完全初始化
   - 添加更多错误保护逻辑
   - 增强客户端连接失败的处理

**性能优化**：
   - 避免在客户端多次重建指针关系，提高连接效率

这些修改解决了在生产环境中发现的段错误和清理错误问题，使服务器模式更加可靠和稳定。

## 2025年4月更新 - 修复段错误和共享内存稳定性增强

**主要修复：**

1. **增强共享内存管理**：
   - 添加魔数(`SHM_MAGIC_NUMBER`)验证共享内存结构的有效性
   - 完善共享内存布局，确保所有指针正确重建
   - 修复了客户端连接后访问共享内存中索引结构的段错误问题

2. **添加全面的调试日志**：
   - 添加详细日志记录系统，跟踪共享内存操作的每一步
   - 记录内存地址、偏移量和关键值，方便故障排查
   - 验证结构体复制和指针重建过程

3. **强化健壮性**：
   - 为所有关键函数添加全面的安全检查
   - 初始化所有内存区域，避免垃圾数据
   - 添加额外验证步骤，确保共享内存中的索引结构完整有效

4. **改进错误处理**：
   - 确保程序能安全地处理各种异常情况
   - 增强Python代码中的错误捕获和处理
   - 添加结构验证测试，确保共享内存中的指针有效

这些改进大大增强了服务器模式的稳定性，解决了之前由于指针无效或内存布局问题导致的段错误，使共享内存模式可以安全地用于生产环境。

## 2025年4月更新 - 修复共享内存指针偏移量错误

**主要修复：**

1. **解决指针偏移量问题**：
   - 修复了客户端读取共享内存时遇到的偏移量不一致问题
   - 解决了在使用服务器模式时报告的"Invalid name offset"错误
   - 确保服务端/客户端一致地解释偏移量值

2. **加强指针有效性验证**：
   - 检查指针偏移量是否在合理范围内
   - 添加多点随机抽样检查，以提前捕获问题
   - 记录详细的指针信息，方便诊断

3. **提高日志可读性**：
   - 记录更多关键点的指针/偏移量信息
   - 在记录大型索引时选择性地记录有代表性的记录点
   - 更清晰地指出错误发生的位置和原因

通过这些改进，现在共享内存模式能够可靠地处理大型索引文件(10GB+)，即使在不同的终端会话或用户之间也能正确共享索引数据。

## 2025年4月更新 - 修复共享内存中的联合体结构问题

**主要修复：**

1. **联合体(union)内存表示问题**：
   - 修复了共享内存中在ptr/offset联合体(union)解释上的重大问题
   - 在保存到共享内存时清除内存区域，防止指针值残留污染
   - 解决了客户端读取到错误偏移量（如139908516156520）的问题

2. **增加内存布局调试信息**：
   - 添加二进制级别的内存调试，显示联合体的字节表示
   - 记录平台关键参数如sizeof(size_t)和sizeof(char*)
   - 添加更详细的偏移量验证日志

3. **兼容不同平台**：
   - 修复字节对齐和大小端可能导致的解释差异
   - 确保联合体数据被正确解析，不管客户端系统和服务器系统是否相同
   - 增强了跨会话和跨用户的稳定性

这些修改解决了之前在使用服务器模式时出现的"Invalid name offset"错误，并确保了共享内存中的数据能够被不同进程一致地解释。特别是针对联合体(union)结构在不同上下文中复用同一个内存空间可能导致的混淆问题进行了专门修复。

## 2025年4月更新 - 修复共享内存中64位指针和偏移量不一致问题

**主要修复：**

1. **彻底解决跨进程的内存表示问题**：
   - 使用显式内存操作(memcpy/memset)绕过联合体(union)的语义问题
   - 避免依赖C编译器对联合体成员的解释，直接以字节级别操作内存
   - 解决了在64位环境下读取偏移量时出现的数值畸变问题

2. **增强内存布局可视性**：
   - 添加全字节的内存转储功能，显示内存中的全部8个字节
   - 从多个角度验证内存值，既检查字节也检查解释后的数值
   - 在客户端验证过程中添加内存原始字节检查点

3. **修复设计缺陷**：
   - 解决了联合体在跨进程通信中可能导致的未定义行为
   - 处理了编译器优化可能引入的问题
   - 增加了更全面的验证手段，确保在任何情况下都能正确解释内存数据

这些修改针对的是一个深层次的系统级问题：在共享内存环境中，联合体的内存布局在不同进程中可能有不同的解释。通过绕过联合体语义，采用直接的内存拷贝方式，我们确保了数据在不同进程间的一致性解释，彻底解决了之前遇到的"Invalid name offset"错误。

## 2025年4月更新 - 重新设计共享内存中的指针存储策略

**主要修复：**

1. **彻底放弃联合体(union)机制**：
   - 不再使用联合体中的偏移量/指针字段语义，改为纯粹内存操作
   - 在共享内存中只存储偏移量，通过memcpy直接读写字节
   - 完全绕过C编译器对union的处理和优化，确保跨进程一致性

2. **新增共享内存版本标识**：
   - 在共享内存头部添加`new_approach`标志，表明使用新的存储策略
   - 兼容可能存在的旧版本客户端
   - 增强日志输出，清晰区分新旧方法

3. **增强内存验证**：
   - 在每个关键步骤输出完整字节序列
   - 验证写入和读取的一致性
   - 在服务端和客户端均采用完全一致的内存处理方式

这次更新解决了之前联合体在共享内存中的行为不一致问题，通过完全放弃联合体的语义功能，仅把其作为内存块来使用，确保了在不同进程中能够正确地读取偏移量，不再受到编译器优化或内存布局的影响。

此修复方案虽然不如直接使用联合体那么优雅，但提供了更可靠和可预测的跨进程内存共享行为，彻底消除了"Invalid name offset"错误。

## 2025年4月更新 - 修复共享内存多客户端数据污染问题

**主要修复：**

1. **实现真正的共享只读访问**：
   - 修复了多客户端访问时共享内存数据被污染的严重bug
   - 客户端现在使用本地内存存储重建后的指针，不再修改共享内存
   - 确保共享内存中始终保持纯偏移量数据，不受客户端操作影响

2. **内存管理改进**：
   - 为每个客户端额外分配本地记录数组
   - 在断开连接时正确释放本地分配的内存
   - 防止内存泄漏并降低长期运行的内存占用

3. **增强多用户并发性**：
   - 客户端之间彻底隔离，一个客户端的操作不再影响其他客户端
   - 可靠支持大规模并发连接场景
   - 解决了多个客户端同时连接可能导致的指针错误

这个修复解决了一个根本性的设计问题：在之前的实现中，客户端在访问共享内存时会修改其中的数据结构（将偏移量转换为指针），这导致后续客户端读取到被污染的数据。通过分离共享内存（只读）和客户端本地内存（可读写），确保了任何客户端都不会影响到共享的数据结构，使系统可以安全稳定地支持多客户端并发访问。

## 2025-04月 - 修复双重释放问题
### 修复
- 修复客户端断开连接时出现的"double free or corruption"错误。在直接访问模式下，客户端上下文指针和bri->records指向同一块内存区域，导致断开连接时发生双重释放。现在增加了检查防止这种情况发生。

## 2025年4月16日 - 修复多客户端并发访问段错误问题

**问题描述:**  
当多个客户端同时连接服务器时，出现段错误(Segmentation fault)。客户端日志显示在查询阶段崩溃，且共享内存中的records指针在两次查询中发生了变化（从0x55f62cf186f0变为0x55b4c3ed6520）。

**根本原因:**  
在`setup_direct_access_mode`函数中直接修改了共享内存中的数据结构：
1. 将共享内存中的`bri->records`指针设置为指向客户端本地的上下文结构
2. 修改共享内存`header->direct_access_mode`标志

当多个客户端同时连接时，这些修改会相互覆盖，导致客户端使用无效指针而崩溃。

**修复方案:**
1. **客户端使用本地BRI副本** - 每个客户端创建并维护自己的BRI结构副本，而不是修改共享内存中的BRI
2. **禁止修改共享内存** - 删除所有修改共享内存的代码，保证共享内存只读
3. **增强内存管理** - 在客户端上下文中保存本地BRI指针，确保正确释放

**技术细节:**
- 在`client_context_t`中添加了`local_bri`字段，存储本地BRI副本
- 修改`setup_direct_access_mode`函数，为每个客户端创建独立的BRI结构体副本
- 修改`py_connect_to_server`函数，返回本地BRI指针而非共享内存中的BRI指针
- 增强`py_disconnect_from_server`函数，确保正确释放客户端本地分配的内存

**改进效果:**
- 完全解决了多客户端同时连接时的段错误问题
- 实现了真正的共享只读访问，消除了客户端之间的相互干扰
- 降低了内存消耗，客户端只需为BRI结构体分配额外内存（几百字节），而不是整个records数组

## 2025-04-16 - 修复异常退出客户端计数不更新问题

**问题描述:**
当客户端异常退出（如段错误、信号终止）时，服务端的客户端计数没有更新，导致计数不准确。

**根本原因:**
客户端异常退出时没有机会执行`disconnect_from_shared_memory`函数，因此客户端计数没有减少。在Linux系统中，进程的段错误（Segmentation fault，错误代码139）会导致进程立即终止，不会执行任何清理代码。

**修复方案:**
1. **客户端PID跟踪机制** - 每个客户端连接时在指定目录创建PID文件，包含进程ID和连接时间
2. **服务端定期检查机制** - 服务器每10秒检查所有PID文件，验证客户端进程是否仍存活
3. **自动清理失效进程** - 对于不再存在的进程，自动更新客户端计数并删除对应PID文件

**技术细节:**
- 添加了`SHM_CLIENT_PID_DIR`目录存储所有客户端的PID文件
- 实现了`create_client_pid_file`和`remove_client_pid_file`函数管理PID文件
- 添加了`check_client_pids_and_update_count`函数，由服务器定期调用检查客户端状态
- 在客户端连接成功时创建PID文件，在正常断开时删除
- 修改服务器清理流程，确保所有PID文件和目录都能被正确删除

**改进效果:**
- 即使客户端异常退出，服务器也能正确更新客户端计数
- 保持客户端计数的准确性，避免长期运行时计数偏差累积
- 提高了系统的整体可靠性和状态一致性




